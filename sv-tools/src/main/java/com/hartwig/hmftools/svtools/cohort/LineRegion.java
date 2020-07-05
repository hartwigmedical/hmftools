package com.hartwig.hmftools.svtools.cohort;

import com.hartwig.hmftools.common.utils.sv.SvRegion;

public class LineRegion
{
    public final SvRegion Region;
    public final LineElementType LineType;
    public int BreakendCount;

    public LineRegion(final SvRegion region, LineElementType type)
    {
        Region = region;
        LineType = type;
        BreakendCount = 1;
    }

    public String toString() { return String.format("loc(%s) type(%s) breakends(%d)", Region, LineType, BreakendCount); }
}
