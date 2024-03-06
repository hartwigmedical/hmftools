package com.hartwig.hmftools.bamtools.common;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public class PartitionTask
{
    public final ChrBaseRegion Region;
    @Nullable
    public final ChrBaseRegion PreviousRegion;
    public final int TaskId;

    public PartitionTask(final ChrBaseRegion region, @Nullable final ChrBaseRegion previousRegion, final int taskId)
    {
        Region = region;
        PreviousRegion = previousRegion;
        TaskId = taskId;
    }
}
