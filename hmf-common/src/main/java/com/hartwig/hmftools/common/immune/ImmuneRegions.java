package com.hartwig.hmftools.common.immune;

import java.sql.Ref;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public final class ImmuneRegions
{
    private static final int TRANSCRIPT_FLANK = 1_000;
    
    public static final ChrBaseRegion IGK_REGION_V37 = new ChrBaseRegion("2", 89156674, 90538397);
    public static final ChrBaseRegion IGH_REGION_V37 = new ChrBaseRegion("14", 106032614, 107288051);
    public static final ChrBaseRegion IGL_REGION_V37 = new ChrBaseRegion("22", 22380474, 23265085);

    public static final ChrBaseRegion IGK_REGION_V38 = new ChrBaseRegion("chr2", 88857161, 90315836);
    public static final ChrBaseRegion IGH_REGION_V38 = new ChrBaseRegion("chr14", 105586437, 106879844);
    public static final ChrBaseRegion IGL_REGION_V38 = new ChrBaseRegion("chr22", 22026076, 22922913);

    public static final ChrBaseRegion TRG_REGION_V37 = new ChrBaseRegion("7", 38279181, 38407483);
    public static final ChrBaseRegion TRB_REGION_V37 = new ChrBaseRegion("7", 141999017, 142510554);
    public static final ChrBaseRegion TRAD_REGION_V37 = new ChrBaseRegion("14", 22090446, 23021099);

    public static final ChrBaseRegion TRG_REGION_V38 = new ChrBaseRegion("chr7", 38239580, 38367882);
    public static final ChrBaseRegion TRB_REGION_V38 = new ChrBaseRegion("chr7", 142299177, 142812869);
    public static final ChrBaseRegion TRAD_REGION_V38 = new ChrBaseRegion("chr14", 21622293, 22552156);
    
    public static final ChrBaseRegion HLA_A_REGION_V37 = new ChrBaseRegion("6", 29909037 - TRANSCRIPT_FLANK, 29913661 + TRANSCRIPT_FLANK);
    public static final ChrBaseRegion HLA_B_REGION_V37 = new ChrBaseRegion("6", 31321649 - TRANSCRIPT_FLANK, 31324964 + TRANSCRIPT_FLANK);
    public static final ChrBaseRegion HLA_C_REGION_V37 = new ChrBaseRegion("6", 31236526 - TRANSCRIPT_FLANK, 31239863 + TRANSCRIPT_FLANK);

    public static final ChrBaseRegion HLA_A_REGION_V38 = new ChrBaseRegion("chr6", 29942532 - TRANSCRIPT_FLANK, 29945870 + TRANSCRIPT_FLANK);
    public static final ChrBaseRegion HLA_B_REGION_V38 = new ChrBaseRegion("chr6", 31353872 - TRANSCRIPT_FLANK, 31367067 + TRANSCRIPT_FLANK);
    public static final ChrBaseRegion HLA_C_REGION_V38 = new ChrBaseRegion("chr6", 31268749 - TRANSCRIPT_FLANK, 31272092 + TRANSCRIPT_FLANK);

    public static final List<ChrBaseRegion> IG_REGIONS_V37 = Lists.newArrayList(IGK_REGION_V37, IGH_REGION_V37, IGL_REGION_V37);
    public static final List<ChrBaseRegion> IG_REGIONS_V38 = Lists.newArrayList(IGK_REGION_V38, IGH_REGION_V38, IGL_REGION_V38);

    public static final List<ChrBaseRegion> TR_REGIONS_V37 = Lists.newArrayList(TRG_REGION_V37, TRB_REGION_V37, TRAD_REGION_V37);
    public static final List<ChrBaseRegion> TR_REGIONS_V38 = Lists.newArrayList(TRG_REGION_V38, TRB_REGION_V38, TRAD_REGION_V38);
    
    public static final List<ChrBaseRegion> HLA_REGIONS_V37 = Lists.newArrayList(HLA_A_REGION_V37, HLA_B_REGION_V37, HLA_C_REGION_V37);
    public static final List<ChrBaseRegion> HLA_REGIONS_V38 = Lists.newArrayList(HLA_A_REGION_V38, HLA_B_REGION_V38, HLA_C_REGION_V38);

    public static List<ChrBaseRegion> getIgRegions(final RefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion.is37() ? IG_REGIONS_V37 : IG_REGIONS_V38;
    }

    public static List<ChrBaseRegion> getTrRegions(final RefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion.is37() ? TR_REGIONS_V37 : TR_REGIONS_V38;
    }
    
    public static List<ChrBaseRegion> getHlaRegions(final RefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion.is37() ? HLA_REGIONS_V37 : HLA_REGIONS_V38;
    }

    public static ChrBaseRegion getIgRegion(final String geneName, final RefGenomeVersion refGenomeVersion)
    {
        if(geneName.equals("IGK"))
            return refGenomeVersion.is37() ? IGK_REGION_V37 : IGK_REGION_V38;
        else if(geneName.equals("IGH"))
            return refGenomeVersion.is37() ? IGH_REGION_V37 : IGH_REGION_V38;
        else if(geneName.equals("IGL"))
            return refGenomeVersion.is37() ? IGL_REGION_V37 : IGL_REGION_V38;
        else
            return null;
    }
    
    public static ChrBaseRegion getHlaRegion(final String geneName, final RefGenomeVersion refGenomeVersion)
    {
        if(geneName.equals("HLA-A"))
            return refGenomeVersion.is37() ? HLA_A_REGION_V37 : HLA_A_REGION_V38;
        else if(geneName.equals("HLA-B"))
            return refGenomeVersion.is37() ? HLA_B_REGION_V37 : HLA_B_REGION_V38;
        else if(geneName.equals("HLA-C"))
            return refGenomeVersion.is37() ? HLA_C_REGION_V37 : HLA_C_REGION_V38;
        else
            return null;
    }
}
