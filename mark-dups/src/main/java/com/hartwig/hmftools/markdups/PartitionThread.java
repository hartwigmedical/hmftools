package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final MarkDupsConfig mConfig;

    private final BamWriter mBamWriter;
    private final Queue<ChrBaseRegion> mPartitions;

    private final PartitionReader mPartitionReader;

    public PartitionThread(
            final int threadId, final MarkDupsConfig config, final Queue<ChrBaseRegion> partitions, final FileWriterCache fileWriterCache,
            final PartitionDataStore partitionDataStore)
    {
        mConfig = config;

        mPartitions = partitions;

        mBamWriter = fileWriterCache.getBamWriter(String.valueOf(threadId));

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
                ChrBaseRegion partition = mPartitions.remove();

                /*
                PartitionReader reader = new PartitionReader(partition.Region, mConfig, mRefSamReader, mNewSamReader, mReadWriter);

                if(partition.TaskId > 0 && (partition.TaskId % 100) == 0)
                {
                    BT_LOGGER.info("processing partition({}), remaining({})", partition.TaskId, mPartitions.size());
                }
                */

                mPartitionReader.setupRegion(partition);
                mPartitionReader.processRegion();

                // reader.run();
                // mStats.merge(reader.stats());
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
