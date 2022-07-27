package com.hartwig.hmftools.common.sv;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public final class ExcludedRegions
{
    public static final String POLY_G_INSERT = "GGGGGGGGGGGGGGGG";
    public static final String POLY_C_INSERT = "CCCCCCCCCCCCCCCC";
    public static final int POLY_G_LENGTH = POLY_G_INSERT.length();

    // regions of high-depth due to issues in the ref-genome or sequencing

    // LINC00486
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_37 = new ChrBaseRegion("2", 33141260, 33141700);
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_38 = new ChrBaseRegion("chr2", 32916190, 32916630);

    public static ChrBaseRegion getPolyGRegion(final RefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion.is37() ? EXCLUDED_REGION_1_REF_37 : EXCLUDED_REGION_1_REF_38;
    }

    public static final List<ChrBaseRegion> POLY_G_REGIONS_V37 = Lists.newArrayList();
    public static final List<ChrBaseRegion> POLY_G_REGIONS_V38 = Lists.newArrayList();

    static
    {
        POLY_G_REGIONS_V37.add(EXCLUDED_REGION_1_REF_37);
        POLY_G_REGIONS_V37.add(new ChrBaseRegion("4", 41218427, 41218467));
        POLY_G_REGIONS_V37.add(new ChrBaseRegion("17", 42646418, 42646458));

        POLY_G_REGIONS_V38.add(EXCLUDED_REGION_1_REF_38);
        POLY_G_REGIONS_V38.add(new ChrBaseRegion("chr4", 41216410, 41216450));
        POLY_G_REGIONS_V38.add(new ChrBaseRegion("chr17", 44569050, 44569090));
    }

    public static List<ChrBaseRegion> getPolyGRegions(final RefGenomeVersion refGenomeVersion)
    {
        return refGenomeVersion.is37() ? POLY_G_REGIONS_V37 : POLY_G_REGIONS_V38;
    }

}
