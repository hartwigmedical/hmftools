package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.svprep.CombinedReadGroups;
import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SvConfig;

public class PartitionThread extends Thread
{
    private final String mChromosome;
    private final SvConfig mConfig;
    private final CombinedReadGroups mCombinedReadGroups;
    private final ResultsWriter mWriter;
    private final CombinedStats mCombinedStats;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final String chromosome, final SvConfig config, final Queue<PartitionTask> partitions,
            final CombinedReadGroups combinedReadGroups, final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mCombinedReadGroups = combinedReadGroups;
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

                PartitionSlicer slicer = new PartitionSlicer(
                        partition.TaskId, partition.Region, mConfig, mCombinedReadGroups, mWriter, mCombinedStats);

                boolean logAndGc = partition.TaskId > 0 && (partition.TaskId % 10) == 0;
                if(logAndGc)
                {
                    SV_LOGGER.info("chromosome({}) processing partition({}), remaining({})",
                            mChromosome, partition.TaskId, mPartitions.size());
                }

                slicer.run();

                if(logAndGc)
                    System.gc();
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
        }
    }
}
