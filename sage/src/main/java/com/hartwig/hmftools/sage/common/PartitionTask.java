package com.hartwig.hmftools.sage.common;

import static java.lang.String.format;

import java.util.Comparator;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class PartitionTask
{
    public final ChrBaseRegion Partition;
    public final int TaskId;

    public PartitionTask(final ChrBaseRegion partition, final int taskId)
    {
        Partition = partition;
        TaskId = taskId;
    }

    public String toString() { return format("id(%d) region(%s)", TaskId, Partition); }

    public static class PartitionTaskComparator implements Comparator<PartitionTask>
    {
        public int compare(final PartitionTask first, final PartitionTask second)
        {
            return first.Partition.compareTo(second.Partition);
        }
    }

}
