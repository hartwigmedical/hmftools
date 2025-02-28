package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FilterConstants
{
    public final int MinQual;
    public final int MinQualHotspot;

    public static final int DEFAULT_MIN_QUAL = 30;
    public static final int DEFAULT_MIN_QUAL_PANEL = 60;

    private static final String CFG_MIN_QUAL = "min_qual";
    private static final String CFG_MIN_QUAL_HOTSPOT = "hotspot_min_qual";

    public final int MinLength;

    public final int MinSupportJunction;
    public final int MinSupportSgl;
    public final int MinSupportHotspot;

    public static final int DEFAULT_MIN_SUPPORT_JUNCTION = 4;
    public static final int DEFAULT_MIN_SUPPORT_SGL = 6;
    public static final int DEFAULT_MIN_SUPPORT_HOTSPOT = 2;

    private static final String CFG_MIN_SUPPORT = "min_support";
    private static final String CFG_MIN_SUPPORT_HOTSPOT = "hotspot_min_support";
    private static final String CFG_MIN_SUPPORT_SGL = "sgl_min_support";

    public final double MinAfJunction;
    public final double MinAfSgl;
    public final double MinAfHotspot;

    public static final double DEFAULT_MIN_AF_SGL = 0.05;
    public static final double DEFAULT_MIN_AF_HOTSPOT = 0.001;
    public static final double DEFAULT_MIN_AF_JUNCTION = 0.002;

    public final boolean FilterSGLs;
    private static final String FILTER_SGLS = "filter_sgls";

    public final List<ChrBaseRegion> PolyGcRegions;
    public final ChrBaseRegion LowQualRegion;

    public final int PonDistance;
    public static final int DEFAULT_PON_DISTANCE = 3;
    public static final String PON_DISTANCE = "pon_distance";

    public static final int SHORT_CALLING_SIZE = 1000;

    public static final int MIN_TRIMMED_ANCHOR_LENGTH = MIN_ANCHOR_LENGTH;
    public static final int MIN_AVG_FRAG_FACTOR = 3;
    public static final double MIN_AVG_FRAG_STD_DEV_FACTOR = 0.6;

    public static final int DEL_ARTEFACT_SHORT_LENGTH = 3000;
    public static final int DEL_ARTEFACT_MAX_HOMOLOGY = 10;
    public static final double DEL_ARTEFACT_LENGTH_FACTOR = 1.5;
    public static final double DEL_ARTEFACT_MIN_AF = 0.05;

    public static final int INV_SHORT_LENGTH = 3000;
    public static final double INV_SHORT_MIN_AF_LOWER = 0.1;
    public static final double INV_SHORT_MIN_AF_HIGHER = 0.2;
    public static final int INV_SHORT_MAX_HOMOLOGY_LOWER = 3;
    public static final int INV_SHORT_MAX_HOMOLOGY_HIGHER = 6;
    public static final int INV_SHORT_RATE_LOWER = 200;
    public static final int INV_SHORT_RATE_HIGHER = 400;

    public static final ChrBaseRegion PMS2_V37 = new ChrBaseRegion("7", 6002870, 6058756); // has 10K buffer
    public static final ChrBaseRegion PMS2_V38 = new ChrBaseRegion("chr7", 5960925, 6019106);

    public static final String PON_INS_SEQ_FWD_STRAND = "GTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA";
    public static final String PON_INS_SEQ_REV_STRAND = "TTTTTAATGATACGGCGACCACCGAGATCTACAC";

    public static final double GERMLINE_AF_THRESHOLD = 0.1;
    public static final double GERMLINE_AD_THRESHOLD = 0.01;

    public static FilterConstants from(final ConfigBuilder configBuilder)
    {
        boolean targetedMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        RefGenomeVersion refGenVersion = RefGenomeVersion.from(configBuilder);

        // use targeted panel defaults where applicable
        boolean filterSgls = targetedMode || configBuilder.hasFlag(FILTER_SGLS);

        int minSupport = configBuilder.getInteger(CFG_MIN_SUPPORT);
        int minSupportHotspot = configBuilder.getInteger(CFG_MIN_SUPPORT_HOTSPOT);
        int minSupportSgl = configBuilder.getInteger(CFG_MIN_SUPPORT_SGL);

        int minQual = configBuilder.hasValue(CFG_MIN_QUAL) || !targetedMode ?
                configBuilder.getInteger(CFG_MIN_QUAL) : DEFAULT_MIN_QUAL_PANEL;

        int minQualHotspot = configBuilder.getInteger(CFG_MIN_QUAL_HOTSPOT);

        int ponDistance = configBuilder.getInteger(PON_DISTANCE);

        return new FilterConstants(
                minQual, minQualHotspot, MIN_VARIANT_LENGTH,
                minSupport, minSupportSgl, minSupportHotspot,
                DEFAULT_MIN_AF_JUNCTION, DEFAULT_MIN_AF_SGL, DEFAULT_MIN_AF_HOTSPOT,
                filterSgls,
                getPolyGRegions(refGenVersion),
                refGenVersion == V37 ? PMS2_V37 : PMS2_V38,
                ponDistance);
    }

    @VisibleForTesting
    public static FilterConstants from(boolean filterSgls, final RefGenomeVersion refGenomeVersion, int ponDistance)
    {
        return new FilterConstants(
                DEFAULT_MIN_QUAL, DEFAULT_MIN_QUAL, MIN_VARIANT_LENGTH,
                DEFAULT_MIN_SUPPORT_JUNCTION, DEFAULT_MIN_SUPPORT_SGL, DEFAULT_MIN_SUPPORT_HOTSPOT,
                DEFAULT_MIN_AF_JUNCTION, DEFAULT_MIN_AF_SGL, DEFAULT_MIN_AF_HOTSPOT,
                filterSgls,
                getPolyGRegions(refGenomeVersion), refGenomeVersion == V37 ? PMS2_V37 : PMS2_V38,
                ponDistance);
    }

    public FilterConstants(
            final int minQual, final int minQualHotspot,
            final int minLength, final int minSupportJunction, final int minSupportSgl, final int minSupportHotspot,
            final double minAfJunction, final double minAfSgl, final double minAfHotspot,
            final boolean filterSGLs, final List<ChrBaseRegion> polyGcRegions, final ChrBaseRegion lowQualRegion, final int ponDistance)
    {
        MinQual = minQual;
        MinQualHotspot = minQualHotspot;
        MinLength = minLength;
        MinSupportJunction = minSupportJunction;
        MinSupportSgl = minSupportSgl;
        MinSupportHotspot = minSupportHotspot;
        MinAfJunction = minAfJunction;
        MinAfSgl = minAfSgl;
        MinAfHotspot = minAfHotspot;
        FilterSGLs = filterSGLs;
        PolyGcRegions = polyGcRegions;
        LowQualRegion = lowQualRegion;
        PonDistance = ponDistance;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(PON_DISTANCE, "PON permitted margin", DEFAULT_PON_DISTANCE);
        configBuilder.addFlag(FILTER_SGLS, "Filter SGLs from VCF, intended for tumor-only mode, default=true in target panel");

        configBuilder.addInteger(CFG_MIN_QUAL,
                "Min qual, panel default: " + String.valueOf(DEFAULT_MIN_QUAL_PANEL), DEFAULT_MIN_QUAL);

        configBuilder.addInteger(CFG_MIN_QUAL_HOTSPOT, "Hotspot min qual", DEFAULT_MIN_QUAL);

        configBuilder.addInteger(CFG_MIN_SUPPORT, "Min support", DEFAULT_MIN_SUPPORT_JUNCTION);
        configBuilder.addInteger(CFG_MIN_SUPPORT_HOTSPOT, "Hotspot min support", DEFAULT_MIN_SUPPORT_HOTSPOT);
        configBuilder.addInteger(CFG_MIN_SUPPORT_SGL, "SGL min support", DEFAULT_MIN_SUPPORT_SGL);
    }
}
