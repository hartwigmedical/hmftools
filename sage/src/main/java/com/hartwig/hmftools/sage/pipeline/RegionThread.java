package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RegionThread extends Thread
{
    private final String mChromosome;
    private final SageConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;

    private final Map<String, QualityRecalibrationMap> mQualityRecalibrationMap;
    private final Coverage mCoverage;
    private  final PhaseSetCounter mPhaseSetCounter;

    private final Queue<PartitionTask> mPartitions;
    private final RegionResults mRegionResults;

    // cache of chromosome-specific ref data
    private final List<BaseRegion> mPanelRegions;
    private final List<VariantHotspot> mHotspots;
    private final List<TranscriptData> mTranscripts;
    private final List<BaseRegion> mHighConfidenceRegions;

    private final Map<String,SamReader> mBamReaders;

    public RegionThread(
            final String chromosome, final SageConfig config,
            final Map<String,QualityRecalibrationMap> qualityRecalibrationMap, final Coverage coverage,
            final PhaseSetCounter phaseSetCounter, final List<BaseRegion> panelRegions, final List<VariantHotspot> hotspots,
            final List<TranscriptData> transcripts, final List<BaseRegion> highConfidenceRegions,
            final Queue<PartitionTask> partitions, final RegionResults regionResults)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mCoverage = coverage;
        mPhaseSetCounter = phaseSetCounter;

        mPanelRegions = panelRegions;
        mHighConfidenceRegions = highConfidenceRegions;
        mHotspots = hotspots;
        mTranscripts = transcripts;

        mRegionResults = regionResults;
        mPartitions = partitions;

        // create readers for each sample and BAM
        mBamReaders = Maps.newHashMap();

        List<String> allSamples = Lists.newArrayList(mConfig.TumorIds);
        allSamples.addAll(mConfig.ReferenceIds);

        List<String> allBams = Lists.newArrayList(mConfig.TumorBams);
        allBams.addAll(mConfig.ReferenceBams);

        for(int i = 0; i < allSamples.size(); i++)
        {
            final String sample = allSamples.get(i);
            final String bamFile = allBams.get(i);

            SamReader bamReader = SamReaderFactory.makeDefault()
                    .validationStringency(mConfig.Stringency)
                    .referenceSource(new ReferenceSource(mRefGenome))
                    .open(new File(bamFile));

            mBamReaders.put(sample, bamReader);
        }

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
        }

        try
        {
            for(SamReader bamReader : mBamReaders.values())
            {
                bamReader.close();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }

    private RegionTask createRegionTask(final PartitionTask partitionTask)
    {
        ChrBaseRegion region = partitionTask.Partition;

        List<BaseRegion> regionPanel = mPanelRegions != null ? mPanelRegions.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                : Lists.newArrayList();

        List<VariantHotspot> regionHotspots = mHotspots != null ? mHotspots.stream()
                .filter(x -> region.containsPosition(x.position())).collect(Collectors.toList()) : Lists.newArrayList();

        List<TranscriptData> regionsTranscripts = mTranscripts != null ? mTranscripts.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.TransStart, x.TransEnd)).collect(Collectors.toList())
                : Lists.newArrayList();

        List<BaseRegion> regionHighConfidence = mHighConfidenceRegions != null ? mHighConfidenceRegions.stream()
                .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                : Lists.newArrayList();

        return new RegionTask(
                partitionTask.TaskId, region, mRegionResults, mConfig, mRefGenome, regionHotspots, regionPanel, regionsTranscripts,
                regionHighConfidence, mQualityRecalibrationMap, mPhaseSetCounter, mCoverage, mBamReaders);
    }
}
