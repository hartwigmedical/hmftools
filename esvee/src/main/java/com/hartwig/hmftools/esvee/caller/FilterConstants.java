package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_FRAG_LOW_VAF;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_ISOLATED;
import static com.hartwig.hmftools.esvee.common.FilterType.INV_SHORT_LOW_VAF_HOM;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.common.FilterType;

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

    public static final int DEFAULT_PON_DISTANCE = 3;
    public static final int DEFAULT_SGL_PON_DISTANCE = 10;

    public static final int PON_SHORT_INDEL_LENGTH = 200;
    public static final int PON_MAX_INS_SEQ_LENGTH = 100;
    public static final int PON_SHORT_INDEL_PON_DISTANCE = 10;
    public static final int PON_SHORT_INDEL_MAX_REPEAT_PON_DISTANCE = 20;
    public static final int PON_SHORT_INDEL_MAX_REPEATS = 11;

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

    public static final int INV_SHORT_FRAGMENT_LENGTH = 300;
    public static final double INV_SHORT_FRAGMENT_MIN_AF = 0.05;
    public static final int INV_SHORT_FRAGMENT_AF_RATIO = 50;

    public static final int PRIME_MAX_PERMITTED_RANGE = 50;
    public static final int PRIME_MAX_BASE_FACTOR = 10;
    public static final int PRIME_MAX_SGL_FACTOR = 5;

    public static final int INV_ADJACENT_LENGTH = 100;
    public static final int INV_ADJACENT_MIN_UPS = 4;
    public static final Set<FilterType> INV_ADJACENT_EXCLUDED_FILTERS = Sets.newHashSet(
            PON, INV_SHORT_FRAG_LOW_VAF, INV_SHORT_LOW_VAF_HOM, INV_SHORT_ISOLATED);

    // SBX-specific filters
    public static final int SBX_STRAND_BIAS_NON_BND_MIN_FRAGS = 10;
    public static final int SBX_HEURISTIC_THRESHOLD = 50;
    public static final int SBX_HEURISTIC_SGL_THRESHOLD = 30;

    public static final int SBX_HEURISTIC_SHORT_LENGTH = 300;
    public static final double SBX_HEURISTIC_ASM_LENGTH_FACTOR = 600;
    public static final int SBX_HEURISTIC_INEXACT_HOM_FACTOR = 2;
    public static final int SBX_HEURISTIC_INEXACT_HOM_MAX = 50;

    public static final int SBX_HEURISTIC_DUP_INS_LENGTH = 40;
    public static final int SBX_HEURISTIC_SHORT_DUP_INS_PENALTY = 20;

    public static final int SBX_HEURISTIC_INV_SHORT_LENGTH = 100;
    public static final int SBX_HEURISTIC_INV_INS_LENGTH = 20;

    public static final int SBX_HEURISTIC_LOCATION_JUNC_PENALTY = 30;
    public static final int SBX_HEURISTIC_BND_LOCATION_JUNC_PENALTY = 20;

    public static final int SBX_HEURISTIC_BND_INS_PENALTY = 20;
    public static final int SBX_HEURISTIC_BND_LINE_BONUS = 20;
    public static final int SBX_HEURISTIC_SGL_LINE_BONUS = 30;

    public static final double SBX_HEURISTIC_BND_LINE_PERC = 0.9;

    public static final ChrBaseRegion PMS2_V37 = new ChrBaseRegion("7", 6002870, 6058756); // has 10K buffer
    public static final ChrBaseRegion PMS2_V38 = new ChrBaseRegion("chr7", 5960925, 6019106);

    public static final ChrBaseRegion PMS2CL_V37 = new ChrBaseRegion("7", 6759759, 6803493);
    public static final ChrBaseRegion PMS2CL_V38 = new ChrBaseRegion("chr7", 6725305, 6761392);

    public static final String PON_INS_SEQ_FWD_STRAND_1 = "GCCGTATCATTAAAAA";
    public static final String PON_INS_SEQ_FWD_STRAND_2 = "GTAGATCTCGGTGGTC";

    public static final String PON_INS_SEQ_REV_STRAND_1 = Nucleotides.reverseComplementBases(PON_INS_SEQ_FWD_STRAND_1);
    public static final String PON_INS_SEQ_REV_STRAND_2 = Nucleotides.reverseComplementBases(PON_INS_SEQ_FWD_STRAND_2);

    public static final int PANEL_INCLUSION_BUFFER = 1000;

    public static FilterConstants from(final ConfigBuilder configBuilder)
    {
        boolean targetedMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        RefGenomeVersion refGenVersion = RefGenomeVersion.from(configBuilder);

        // use targeted panel defaults where applicable
        boolean filterSgls = configBuilder.hasFlag(FILTER_SGLS);

        int minSupport = configBuilder.getInteger(CFG_MIN_SUPPORT);
        int minSupportHotspot = configBuilder.getInteger(CFG_MIN_SUPPORT_HOTSPOT);
        int minSupportSgl = configBuilder.getInteger(CFG_MIN_SUPPORT_SGL);

        int minQual = configBuilder.hasValue(CFG_MIN_QUAL) || !targetedMode ?
                configBuilder.getInteger(CFG_MIN_QUAL) : DEFAULT_MIN_QUAL_PANEL;

        int minQualHotspot = configBuilder.getInteger(CFG_MIN_QUAL_HOTSPOT);

        return new FilterConstants(
                minQual, minQualHotspot, MIN_VARIANT_LENGTH,
                minSupport, minSupportSgl, minSupportHotspot,
                DEFAULT_MIN_AF_JUNCTION, DEFAULT_MIN_AF_SGL, DEFAULT_MIN_AF_HOTSPOT,
                filterSgls,
                getPolyGRegions(refGenVersion),
                refGenVersion == V37 ? PMS2_V37 : PMS2_V38);
    }

    @VisibleForTesting
    public static FilterConstants from(boolean filterSgls, final RefGenomeVersion refGenomeVersion, int ponDistance)
    {
        return new FilterConstants(
                DEFAULT_MIN_QUAL, DEFAULT_MIN_QUAL, MIN_VARIANT_LENGTH,
                DEFAULT_MIN_SUPPORT_JUNCTION, DEFAULT_MIN_SUPPORT_SGL, DEFAULT_MIN_SUPPORT_HOTSPOT,
                DEFAULT_MIN_AF_JUNCTION, DEFAULT_MIN_AF_SGL, DEFAULT_MIN_AF_HOTSPOT,
                filterSgls,
                getPolyGRegions(refGenomeVersion), refGenomeVersion == V37 ? PMS2_V37 : PMS2_V38);
    }

    public FilterConstants(
            final int minQual, final int minQualHotspot,
            final int minLength, final int minSupportJunction, final int minSupportSgl, final int minSupportHotspot,
            final double minAfJunction, final double minAfSgl, final double minAfHotspot,
            final boolean filterSGLs, final List<ChrBaseRegion> polyGcRegions, final ChrBaseRegion lowQualRegion)
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
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(FILTER_SGLS, "Filter SGLs from VCF, intended for tumor-only mode, default=true in target panel");

        configBuilder.addInteger(CFG_MIN_QUAL,
                "Min qual, panel default: " + String.valueOf(DEFAULT_MIN_QUAL_PANEL), DEFAULT_MIN_QUAL);

        configBuilder.addInteger(CFG_MIN_QUAL_HOTSPOT, "Hotspot min qual", DEFAULT_MIN_QUAL);

        configBuilder.addInteger(CFG_MIN_SUPPORT, "Min support", DEFAULT_MIN_SUPPORT_JUNCTION);
        configBuilder.addInteger(CFG_MIN_SUPPORT_HOTSPOT, "Hotspot min support", DEFAULT_MIN_SUPPORT_HOTSPOT);
        configBuilder.addInteger(CFG_MIN_SUPPORT_SGL, "SGL min support", DEFAULT_MIN_SUPPORT_SGL);
    }
}
