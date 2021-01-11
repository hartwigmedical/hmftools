package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class LineRegion
{
    public final BaseRegion Region;
    public final LineElementType LineType;
    public int BreakendCount;

    public LineRegion(final BaseRegion region, LineElementType type)
    {
        Region = region;
        LineType = type;
        BreakendCount = 1;
    }

    public String toString() { return String.format("loc(%s) type(%s) breakends(%d)", Region, LineType, BreakendCount); }
}
