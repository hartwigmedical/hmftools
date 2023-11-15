package com.hartwig.hmftools.gripss.pon;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PonSvRegion implements Comparable<PonSvRegion>
{
    public final BaseRegion RegionStart;
    public final Byte OrientStart;
    public final ChrBaseRegion RegionEnd;
    public final Byte OrientEnd;
    public final int PonCount;

    public PonSvRegion(
            final BaseRegion regionStart, final Byte orientStart, final ChrBaseRegion regionEnd, final Byte orientEnd, final int ponCount)
    {
        RegionStart = regionStart;
        OrientStart = orientStart;
        RegionEnd = regionEnd;
        OrientEnd = orientEnd;
        PonCount = ponCount;
    }

    public boolean matches(final BaseRegion svStart, final ChrBaseRegion svEnd, byte orientStart, byte orientEnd)
    {
        return RegionStart.overlaps(svStart) && RegionEnd.overlaps(svEnd) && OrientStart == orientStart && OrientEnd == orientEnd;
    }

    @Override
    public int compareTo(final PonSvRegion other)
    {
        if(RegionStart.start() == other.RegionStart.start())
        {
            if(RegionStart.end() == other.RegionStart.end())
                return 0;

            return RegionStart.end() < other.RegionStart.end() ? -1 : 1;
        }

        return RegionStart.start() < other.RegionStart.start() ? -1 : 1;
    }

    public String toString()
    {
        return String.format("region(%s - %s) orients(%d - %d) pon(%d)", RegionStart, RegionEnd, OrientStart, OrientEnd, PonCount);
    }
}
