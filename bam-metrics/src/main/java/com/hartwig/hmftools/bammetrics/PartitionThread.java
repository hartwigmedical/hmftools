package com.hartwig.hmftools.bammetrics;

import static com.hartwig.hmftools.bammetrics.BmConfig.BM_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.samtools.BamSlicer;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final String mChromosome;
    private final BmConfig mConfig;
    private final Metrics mCombinedMetrics;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final String chromosome, final BmConfig config, final Queue<PartitionTask> partitions, final Metrics combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mCombinedMetrics = combinedStats;
        mPartitions = partitions;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
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

                PartitionSlicer slicer = new PartitionSlicer(
                        partition.Region, mConfig, mSamReader, mBamSlicer, mCombinedMetrics);

                boolean logAndGc = partition.TaskId > 0 && (partition.TaskId % 10) == 0;

                if(logAndGc)
                {
                    BM_LOGGER.debug("chromosome({}) processing partition({}), remaining({})",
                            mChromosome, partition.TaskId, mPartitions.size());
                }

                slicer.run();

                if(logAndGc)
                    System.gc();
            }
            catch(NoSuchElementException e)
            {
                BM_LOGGER.trace("all tasks complete");
                break;
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            BM_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
