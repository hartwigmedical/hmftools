package com.hartwig.hmftools.sage.pipeline;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ChromosomePipeline implements AutoCloseable
{
    private final String mChromosome;
    private final SageCallConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;

    private final Map<String, BqrRecordMap> mQualityRecalibrationMap;
    private final MsiJitterCalcs mMsiJitterCalcs;
    private final Coverage mCoverage;
    private  final PhaseSetCounter mPhaseSetCounter;

    private final VcfWriter mVcfWriter;
    private final FragmentLengths mFragmentLengths;
    private final Queue<PartitionTask> mPartitions;
    private final RegionResults mRegionResults;

    // cache of chromosome-specific ref data
    private final List<BaseRegion> mPanelRegions;
    private final List<SimpleVariant> mHotspots;
    private final List<TranscriptData> mTranscripts;
    private final List<BaseRegion> mHighConfidenceRegions;

    public ChromosomePipeline(
            final String chromosome, final SageCallConfig config,
            final ReferenceData refData, final Map<String, BqrRecordMap> qualityRecalibrationMap, final MsiJitterCalcs msiJitterCalcs,
            final Coverage coverage, final PhaseSetCounter phaseSetCounter, final VcfWriter vcfWriter, final FragmentLengths fragmentLengths)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.Common.RefGenomeFile);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mMsiJitterCalcs = msiJitterCalcs;
        mCoverage = coverage;
        mPhaseSetCounter = phaseSetCounter;

        mVcfWriter = vcfWriter;
        mFragmentLengths = fragmentLengths;

        final Chromosome chr = HumanChromosome.contains(chromosome)
                ? HumanChromosome.fromString(chromosome) : MitochondrialChromosome.fromString(chromosome);

        mPanelRegions = refData.PanelWithHotspots.get(chr);
        mHotspots = refData.Hotspots.get(chr);
        mTranscripts = refData.ChromosomeTranscripts.get(chromosome);
        mHighConfidenceRegions = refData.HighConfidence.get(chr);

        mPartitions = new ConcurrentLinkedQueue<>();
        mRegionResults = new RegionResults(vcfWriter);

        // split chromosome into partitions, filtering for the panel if in use
        ChromosomePartition chrPartition = new ChromosomePartition(config.Common, mRefGenome);
        List<ChrBaseRegion> partitionedRegions = chrPartition.partition(mChromosome);

        int taskId = 0;
        for(int i = 0; i < partitionedRegions.size(); ++i)
        {
            ChrBaseRegion region = partitionedRegions.get(i);

            List<BaseRegion> regionPanel = mPanelRegions != null ? mPanelRegions.stream()
                    .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                    : Lists.newArrayList();

            if(mConfig.PanelOnly && regionPanel.isEmpty())
                continue;

            mPartitions.add(new PartitionTask(region, taskId++));
        }
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public void process()
    {
        int regionCount = mPartitions.size();
        SG_LOGGER.info("chromosome({}) executing {} regions", mChromosome, regionCount);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(mPartitions.size(), mConfig.Common.Threads); ++i)
        {
            workers.add(new RegionThread(
                    mChromosome, mConfig, mQualityRecalibrationMap, mMsiJitterCalcs, mCoverage, mPhaseSetCounter,
                    mPanelRegions, mHotspots, mTranscripts, mHighConfidenceRegions, mPartitions, mRegionResults, mFragmentLengths));
        }

        if(!runThreadTasks(workers))
            System.exit(1);

        SG_LOGGER.debug("chromosome({}) {} regions complete, initial candidates({}) final variants({}) reads({})",
                mChromosome, regionCount, mRegionResults.totalCandidates(), mRegionResults.totalVariants(), mRegionResults.totalReads());

        mVcfWriter.flushChromosome();

        if(mConfig.Common.logPerfStats())
        {
            mRegionResults.logPerfCounters();
            SG_LOGGER.debug("chromosome({}) evidence stats: {}", mChromosome, mRegionResults.evidenceStats().toString());
        }

        if(mConfig.Common.SyncFragments)
            mRegionResults.logSynCounts();

        SG_LOGGER.info("chromosome({}) analysis complete", mChromosome);
    }

    @Override
    public void close() throws IOException
    {
        mRefGenome.close();
    }

}
