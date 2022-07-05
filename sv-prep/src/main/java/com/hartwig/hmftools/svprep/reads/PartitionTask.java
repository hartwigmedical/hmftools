package com.hartwig.hmftools.svprep.reads;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class PartitionTask
{
    public final ChrBaseRegion Region;
    public final int TaskId;

    public PartitionTask(final ChrBaseRegion region, final int taskId)
    {
        Region = region;
        TaskId = taskId;
    }
}
