package com.hartwig.hmftools.esvee.caller.annotation;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;

public class PonSglRegion implements Comparable<PonSglRegion>
{
    public final BaseRegion Region;
    public final Orientation Orient;
    public final int PonCount;

    public PonSglRegion(final BaseRegion region, final Orientation orient, final int ponCount)
    {
        Region = region;
        Orient = orient;
        PonCount = ponCount;
    }

    public boolean matches(final BaseRegion svRegion, Orientation orientation)
    {
        return Region.overlaps(svRegion) && Orient == orientation;
    }

    public String toString()
    {
        return String.format("region(%s) orient(%d) pon(%d)", Region, Orient.asByte(), PonCount);
    }

    @Override
    public int compareTo(final PonSglRegion other)
    {
        if(Region.start() == other.Region.start())
        {
            if(Region.end() == other.Region.end())
                return 0;

            return Region.end() < other.Region.end() ? -1 : 1;
        }

        return Region.start() < other.Region.start() ? -1 : 1;
    }

}
