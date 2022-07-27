package com.hartwig.hmftools.svprep.tools;

import static java.lang.String.format;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class HighDepthRegion
{
    public final ChrBaseRegion Region;
    public int Depth;

    public HighDepthRegion(final ChrBaseRegion region)
    {
        Region = region;
        Depth = 0;
    }

    public String toString() { return format("region(%s) count(%d)", Region, Depth); }
}
