package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.markdups.write.BamWriter;
import com.hartwig.hmftools.markdups.write.FileWriterCache;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final MarkDupsConfig mConfig;

    private final BamWriter mBamWriter;
    private final Queue<ChrBaseRegion> mPartitions;
    private final int mPartitionCount;

    private final PartitionReader mPartitionReader;

    public PartitionThread(
            final int threadId, final MarkDupsConfig config, final Queue<ChrBaseRegion> partitions, final FileWriterCache fileWriterCache,
            final PartitionDataStore partitionDataStore)
    {
        mConfig = config;

        mPartitions = partitions;
        mPartitionCount = partitions.size();

        mBamWriter = fileWriterCache.getPartitionBamWriter(String.valueOf(threadId));

        SamReader samReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mPartitionReader = new PartitionReader(config, samReader, mBamWriter, partitionDataStore);

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
                    MD_LOGGER.info("processed {} partitions, remaining({})", processedCount, remainingCount);
                }

                mPartitionReader.setupRegion(partition);
                mPartitionReader.processRegion();
            }
            catch(NoSuchElementException e)
            {
                MD_LOGGER.trace("all tasks complete");
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
