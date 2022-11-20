package com.hartwig.hmftools.bamtools;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.bamtools.metrics.CombinedStats;
import com.hartwig.hmftools.bamtools.slice.SliceWriter;
import com.hartwig.hmftools.common.samtools.BamSlicer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final BmConfig mConfig;
    private final CombinedStats mCombinedStats;
    private final SliceWriter mSliceWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final BmConfig config, final Queue<PartitionTask> partitions, final CombinedStats combinedStats, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mCombinedStats = combinedStats;
        mSliceWriter = sliceWriter;
        mPartitions = partitions;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();

                PartitionSlicer slicer = new PartitionSlicer(partition.Region, mConfig, mSamReader, mBamSlicer, mCombinedStats, mSliceWriter);

                boolean logAndGc = partition.TaskId > 0 && (partition.TaskId % 10) == 0;

                if(logAndGc)
                {
                    BmConfig.BM_LOGGER.debug("processing partition({}), remaining({})", partition.TaskId, mPartitions.size());
                }

                slicer.run();

                if(logAndGc)
                    System.gc();
            }
            catch(NoSuchElementException e)
            {
                BmConfig.BM_LOGGER.trace("all tasks complete");
                break;
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            BmConfig.BM_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
