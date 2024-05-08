package com.hartwig.hmftools.esvee.caller;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class KnownHotspot
{
    public final ChrBaseRegion RegionStart;
    public final Orientation OrientStart;
    public final ChrBaseRegion RegionEnd;
    public final Orientation OrientEnd;
    public final String GeneInfo;

    public KnownHotspot(
            final ChrBaseRegion regionStart, final Orientation orientStart, final ChrBaseRegion regionEnd, final Orientation orientEnd, final String geneInfo)
    {
        RegionStart = regionStart;
        OrientStart = orientStart;
        RegionEnd = regionEnd;
        OrientEnd = orientEnd;
        GeneInfo = geneInfo;
    }

    public boolean matches(final String chrStart, final String chrEnd, int posStart, int posEnd, Orientation orientStart, Orientation orientEnd)
    {
        if(RegionStart.containsPosition(chrStart, posStart) && OrientStart == orientStart
        && RegionEnd.containsPosition(chrEnd, posEnd) && OrientEnd == orientEnd)
        {
            return true;
        }

        if(RegionStart.containsPosition(chrEnd, posEnd) && OrientStart == orientEnd
        && RegionEnd.containsPosition(chrStart, posStart) && OrientEnd == orientStart)
        {
            return true;
        }

        return false;
    }

}
