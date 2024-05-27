package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.FileWriterCache;

public class PartitionThread extends Thread
{
    private final ReduxConfig mConfig;

    private final BamReader mBamReader;
    private final BamWriter mBamWriter;
    private final Queue<ChrBaseRegion> mPartitions;
    private final int mPartitionCount;

    private final PartitionReader mPartitionReader;

    public PartitionThread(
            final int threadId, final ReduxConfig config, final Queue<ChrBaseRegion> partitions, final FileWriterCache fileWriterCache,
            final PartitionDataStore partitionDataStore)
    {
        mConfig = config;

        mPartitions = partitions;
        mPartitionCount = partitions.size();

        mBamWriter = fileWriterCache.getPartitionBamWriter(String.valueOf(threadId));

        mBamReader = new BamReader(config);

        mPartitionReader = new PartitionReader(config, mBamReader, mBamWriter, partitionDataStore);

        start();
    }

    public PartitionReader partitionReader() { return mPartitionReader; }
    public BamWriter bamWriter() { return mBamWriter; }

    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mPartitions.size();
                int processedCount = mPartitionCount - remainingCount;

                ChrBaseRegion partition = mPartitions.remove();

                if(processedCount > 0 && (processedCount % 100) == 0)
                {
                    RD_LOGGER.info("processed {} partitions, remaining({})", processedCount, remainingCount);
                }

                mPartitionReader.setupRegion(partition);
                mPartitionReader.processRegion();
            }
            catch(NoSuchElementException e)
            {
                RD_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
}
