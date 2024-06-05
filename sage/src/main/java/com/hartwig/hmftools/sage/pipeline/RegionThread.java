package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RegionThread extends Thread
{
    private final String mChromosome;
    private final SageCallConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenomeFile;
    private final RefGenomeSource mRefGenome;

    private final Map<String, BqrRecordMap> mQualityRecalibrationMap;
    private final MsiJitterCalcs mMsiJitterCalcs;
    private final Coverage mCoverage;
    private  final PhaseSetCounter mPhaseSetCounter;

    private final Queue<PartitionTask> mPartitions;
    private final RegionResults mRegionResults;

    // cache of chromosome-specific ref data
    private final List<BaseRegion> mPanelRegions;
    private final List<SimpleVariant> mHotspots;
    private final List<TranscriptData> mTranscripts;
    private final List<BaseRegion> mHighConfidenceRegions;

    private final SamSlicerFactory mSamSlicerFactory;
    private final FragmentLengths mFragmentLengths;

    public RegionThread(
            final String chromosome, final SageCallConfig config,
            final Map<String, BqrRecordMap> qualityRecalibrationMap, final MsiJitterCalcs msiJitterCalcs, final Coverage coverage,
            final PhaseSetCounter phaseSetCounter, final List<BaseRegion> panelRegions, final List<SimpleVariant> hotspots,
            final List<TranscriptData> transcripts, final List<BaseRegion> highConfidenceRegions,
            final Queue<PartitionTask> partitions, final RegionResults regionResults, final FragmentLengths fragmentLengths)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSamSlicerFactory = new SamSlicerFactory();
        mRefGenomeFile = loadRefGenome(config.Common.RefGenomeFile);
        mRefGenome = new RefGenomeSource(mRefGenomeFile);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mMsiJitterCalcs = msiJitterCalcs;
        mCoverage = coverage;
        mPhaseSetCounter = phaseSetCounter;
        mFragmentLengths = fragmentLengths;

        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;
        mHotspots = hotspots;
        mTranscripts = transcripts;

        mRegionResults = regionResults;
        mPartitions = partitions;

        // create readers for each sample and BAM
        mSamSlicerFactory.buildBamReaders(mConfig.TumorIds, mConfig.TumorBams, mConfig.Common, mRefGenomeFile);

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();
                RegionTask task = createRegionTask(partition);

                if(partition.TaskId > 0 && (partition.TaskId % 100) == 0)
                {
                    SG_LOGGER.debug("chromosome({}) regions assigned({}) remaining({})",
                            mChromosome, partition.TaskId, mPartitions.size());
                }

                task.run();
            }
            catch(NoSuchElementException e)
            {
                SG_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                SG_LOGGER.error("thread execution error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        mSamSlicerFactory.close();
    }

    private RegionTask createRegionTask(final PartitionTask partitionTask)
    {
        ChrBaseRegion region = partitionTask.Partition;

        List<BaseRegion> regionPanel = mPanelRegions != null ? mPanelRegions.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                : Lists.newArrayList();

        List<SimpleVariant> regionHotspots = mHotspots != null ? mHotspots.stream()
                .filter(x -> region.containsPosition(x.position())).collect(Collectors.toList()) : Lists.newArrayList();

        List<TranscriptData> regionsTranscripts = mTranscripts != null ? mTranscripts.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.TransStart, x.TransEnd)).collect(Collectors.toList())
                : Lists.newArrayList();

        List<BaseRegion> regionHighConfidence = mHighConfidenceRegions != null ? mHighConfidenceRegions.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                : Lists.newArrayList();

        return new RegionTask(
                partitionTask.TaskId, region, mRegionResults, mConfig, mRefGenome, regionHotspots, regionPanel, regionsTranscripts,
                regionHighConfidence, mQualityRecalibrationMap, mMsiJitterCalcs, mPhaseSetCounter, mCoverage, mSamSlicerFactory,
                mFragmentLengths);
    }
}
