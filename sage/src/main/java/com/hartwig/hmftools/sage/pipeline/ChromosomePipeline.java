package com.hartwig.hmftools.sage.pipeline;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ChromosomePipeline implements AutoCloseable
{
    private final String mChromosome;
    private final SageConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;

    private final Map<String,QualityRecalibrationMap> mQualityRecalibrationMap;
    private final Coverage mCoverage;
    private  final PhaseSetCounter mPhaseSetCounter;

    private final VcfWriter mVcfWriter;
    private final Queue<PartitionTask> mPartitions;
    private final RegionResults mRegionResults;

    // cache of chromosome-specific ref data
    private final List<BaseRegion> mPanelRegions;
    private final List<VariantHotspot> mHotspots;
    private final List<TranscriptData> mTranscripts;
    private final List<BaseRegion> mHighConfidenceRegions;

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    public ChromosomePipeline(
            final String chromosome, final SageConfig config,
            final ReferenceData refData, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap,
            final Coverage coverage, final PhaseSetCounter phaseSetCounter, final VcfWriter vcfWriter)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mCoverage = coverage;
        mPhaseSetCounter = phaseSetCounter;

        mVcfWriter = vcfWriter;

        final Chromosome chr = HumanChromosome.contains(chromosome)
                ? HumanChromosome.fromString(chromosome) : MitochondrialChromosome.fromString(chromosome);

        mPanelRegions = refData.PanelWithHotspots.get(chr);
        mHotspots = refData.Hotspots.get(chr);
        mTranscripts = refData.ChromosomeTranscripts.get(chromosome);
        mHighConfidenceRegions = refData.HighConfidence.get(chr);

        mPartitions = new ConcurrentLinkedQueue<>();
        mRegionResults = new RegionResults(vcfWriter);

        // split chromosome into partitions, filtering for the panel if in use
        ChromosomePartition chrPartition = new ChromosomePartition(config, mRefGenome);
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

        for(int i = 0; i < min(mPartitions.size(), mConfig.Threads); ++i)
        {
            workers.add(new RegionThread(
                    mChromosome, mConfig, mQualityRecalibrationMap, mCoverage, mPhaseSetCounter,
                    mPanelRegions, mHotspots, mTranscripts, mHighConfidenceRegions,
                    mPartitions, mRegionResults));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                SG_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        SG_LOGGER.debug("chromosome({}) {} regions complete, processed {} reads, writing {} variants",
                mChromosome, regionCount, mRegionResults.totalReads(), mRegionResults.totalVariants());

        mVcfWriter.flushChromosome();

        if(mConfig.logPerfStats())
        {
            mRegionResults.logPerfCounters();
            SG_LOGGER.debug("chromosome({}) max memory({})", mChromosome, mRegionResults.maxMemoryUsage());
        }

        SG_LOGGER.info("chromosome({}) analysis complete", mChromosome);
    }

    public int maxMemoryUsage() { return mRegionResults.maxMemoryUsage(); }

    @Override
    public void close() throws IOException
    {
        mRefGenome.close();
    }

}
