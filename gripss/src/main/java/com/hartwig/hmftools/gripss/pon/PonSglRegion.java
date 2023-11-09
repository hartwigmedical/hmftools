package com.hartwig.hmftools.gripss.pon;

import com.hartwig.hmftools.common.region.BaseRegion;

public class PonSglRegion
{
    public final BaseRegion Region;
    public final Byte Orient;
    public final int PonCount;

    public PonSglRegion(final BaseRegion region, final Byte orient, final int ponCount)
    {
        Region = region;
        Orient = orient;
        PonCount = ponCount;
    }

    public boolean matches(final BaseRegion svRegion, byte orientation)
    {
        return Region.overlaps(svRegion) && Orient == orientation;
    }

    public String toString()
    {
        return String.format("region(%s) orient(%d) pon(%d)", Region, Orient, PonCount);
    }

}
