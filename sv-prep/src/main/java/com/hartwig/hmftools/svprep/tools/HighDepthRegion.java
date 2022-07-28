package com.hartwig.hmftools.svprep.tools;

import static java.lang.String.format;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class HighDepthRegion
{
    public final ChrBaseRegion Region;
    public int DepthMin;
    public int DepthMax;

    public HighDepthRegion(final ChrBaseRegion region)
    {
        Region = region;
        DepthMin = 0;
        DepthMax = 0;
    }

    public String toString() { return format("region(%s) depth(min=%d max=%d)", Region, DepthMin, DepthMax); }
}
