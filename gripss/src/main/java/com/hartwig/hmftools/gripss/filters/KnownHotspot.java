package com.hartwig.hmftools.gripss.filters;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class KnownHotspot
{
    public final ChrBaseRegion RegionStart;
    public final Byte OrientStart;
    public final ChrBaseRegion RegionEnd;
    public final Byte OrientEnd;
    public final String GeneInfo;

    public KnownHotspot(
            final ChrBaseRegion regionStart, final Byte orientStart, final ChrBaseRegion regionEnd, final Byte orientEnd, final String geneInfo)
    {
        RegionStart = regionStart;
        OrientStart = orientStart;
        RegionEnd = regionEnd;
        OrientEnd = orientEnd;
        GeneInfo = geneInfo;
    }

    public boolean matches(final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd)
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
