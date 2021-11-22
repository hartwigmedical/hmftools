package com.hartwig.hmftools.gripss.filters;

import com.hartwig.hmftools.common.sv.StructuralVariant;
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

    public boolean matches(final StructuralVariant sv)
    {
        return RegionStart.containsPosition(sv.chromosome(true), sv.position(true).intValue())
            && OrientStart == sv.orientation(true)
            && RegionEnd.containsPosition(sv.chromosome(false), sv.position(false).intValue())
            && OrientEnd == sv.orientation(false);
    }
}
