package com.hartwig.hmftools.bamtools.compare;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import com.hartwig.hmftools.bamtools.common.PartitionTask;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PartitionThread extends Thread
{
    private final CompareConfig mConfig;

    private final SamReader mRefSamReader;
    private final SamReader mNewSamReader;
    private final ReadWriter mReadWriter;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final CompareConfig config, final Queue<PartitionTask> partitions, final ReadWriter readWriter)
    {
        mConfig = config;
        mReadWriter = readWriter;
        mPartitions = partitions;

        mRefSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.RefBamFile));

        mNewSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.NewBamFile));

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();

                PartitionReader reader = new PartitionReader(partition.Region, mConfig, mRefSamReader, mNewSamReader, mReadWriter);

                boolean logAndGc = partition.TaskId > 0 && (partition.TaskId % 10) == 0;

                if(logAndGc)
                {
                    BT_LOGGER.info("processing partition({}), remaining({})", partition.TaskId, mPartitions.size());
                }

                reader.run();

                if(logAndGc)
                    System.gc();
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
            mRefSamReader.close();
            mNewSamReader.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
