package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final MetricsConfig mConfig;
    private final CombinedStats mCombinedStats;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final Queue<PartitionTask> mPartitions;
    private final MetricsWriter mMetricsWriter;

    public PartitionThread(
            final MetricsConfig config, final Queue<PartitionTask> partitions, final CombinedStats combinedStats,
            final MetricsWriter metricsWriter)
    {
        mConfig = config;
        mCombinedStats = combinedStats;
        mPartitions = partitions;
        mMetricsWriter = metricsWriter;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepUnmapped();
    }

    @Override
    public void run()
    {
        if(mPartitions.isEmpty())
            return;

        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();

                BamReader bamReader = new BamReader(partition.Region, mConfig, mSamReader, mBamSlicer, mCombinedStats);

                if(partition.TaskId > 0 && (partition.TaskId % 10) == 0)
                {
                    BT_LOGGER.info("processing partition({}), remaining({})", partition.TaskId, mPartitions.size());
                }

                bamReader.run();

                TargetRegionStats.writeStatistics(mMetricsWriter.targetRegionsWriter(), bamReader.targetRegionStats());

                if(mMetricsWriter.offTargetHighFragmentOverlapWriter() != null)
                {
                    OffTargetFragments.writeEnrichedRegions(
                            mMetricsWriter.offTargetHighFragmentOverlapWriter(), bamReader.offTargetFragments().enrichedFragmentSites());
                }
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
