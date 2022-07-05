package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SvConfig;

public class PartitionThread extends Thread
{
    private final String mChromosome;
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;
    private final CombinedStats mCombinedStats;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final String chromosome, final SvConfig config, final Queue<PartitionTask> partitions, final ResultsWriter writer,
            final CombinedStats combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mWriter = writer;
        mCombinedStats = combinedStats;
        mPartitions = partitions;

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();

                PartitionSlicer slicer = new PartitionSlicer(partition.TaskId, partition.Region, mConfig, mWriter, mCombinedStats);

                if(partition.TaskId > 0 && (partition.TaskId % 100) == 0)
                {
                    SV_LOGGER.debug("chromosome({}) partitions processed({}) remaining({})",
                            mChromosome, partition.TaskId, mPartitions.size());
                }

                slicer.run();
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
        }

        // mSamSlicerFactory.close();
    }

}
