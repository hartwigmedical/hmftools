package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class LineRegion
{
    public final ChrBaseRegion Region;
    public final LineElementType LineType;
    public int BreakendCount;

    public LineRegion(final ChrBaseRegion region, LineElementType type)
    {
        Region = region;
        LineType = type;
        BreakendCount = 1;
    }

    public String toString() { return String.format("loc(%s) type(%s) breakends(%d)", Region, LineType, BreakendCount); }
}
