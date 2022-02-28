package com.hartwig.hmftools.sage.common;

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
}
