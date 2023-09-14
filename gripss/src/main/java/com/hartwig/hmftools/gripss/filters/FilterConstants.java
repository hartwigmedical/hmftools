package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegions;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class FilterConstants
{
    public final int MinTumorQual;
    public final int HardMaxNormalAbsoluteSupport;
    public final double HardMaxNormalRelativeSupport;
    public final double SoftMaxNormalRelativeSupport;
    public final double MinNormalCoverage;
    public final double MinTumorAfBreakend;
    public final double MinTumorAfBreakpoint;
    public final double MaxShortStrandBias;
    public final int MinQualBreakend;
    public final int MinQualBreakpoint;
    public final int MinQualRescueLine;
    public final int MaxHomLengthShortInv;
    public final int MinLength;
    public final int PonDistance;
    public final List<ChrBaseRegion> PolyGcRegions;
    public final ChrBaseRegion LowQualRegion;
    public final boolean FilterSGLs;

    // filters which only apply when reference is present:
    // minNormalCoverage, minRelativeCoverage, maxNormalSupport, shortSRNormalSupport, discordantPairSupport

    // default filter values
    public static final int SHORT_CALLING_SIZE = 1000;
    public static final int HOM_INV_LENGTH = 50;

    public static final int SGL_INS_SEQ_MIN_LENGTH = 16;
    public static final double SGL_MIN_STRAND_BIAS = 0.05;
    public static final double SGL_MAX_STRAND_BIAS = 0.95;

    public static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 100;
    public static final int DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT = 3;
    public static final double DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT = 0.08;
    public static final double DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT = 0.03;
    public static final double DEFAULT_MIN_NORMAL_COVERAGE = 8;
    public static final double DEFAULT_MIN_TUMOR_AF = 0.005;
    public static final double DEFAULT_MIN_TUMOR_AF_SGL = 0.015;
    public static final double DEFAULT_MAX_SHORT_STRAND_BIAS = 0.95;
    public static final int DEFAULT_MIN_QUAL_BREAK_END = 500;
    public static final int DEFAULT_MIN_QUAL_BREAK_POINT = 400;
    public static final int DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION = 500;
    public static final int DEFAULT_MAX_HOM_LENGTH_SHORT_INV = 6;
    public static final int DEFAULT_MIN_LENGTH = 32;
    public static final int DEFAULT_PON_DISTANCE = 3;

    public static final ChrBaseRegion PMS2_V37 = new ChrBaseRegion("7", 6002870, 6058756); // has 10K buffer
    public static final ChrBaseRegion PMS2_V38 = new ChrBaseRegion("chr7", 5960925, 6019106);

    // config overrides
    public static final String HARD_MIN_TUMOR_QUAL_CFG = "hard_min_tumor_qual";
    public static final String HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_CFG = "hard_max_normal_absolute_support";
    public static final String HARD_MAX_NORMAL_RELATIVE_SUPPORT_CFG = "hard_max_normal_relative_support";
    public static final String SOFT_MAX_NORMAL_RELATIVE_SUPPORT_CFG = "soft_max_normal_relative_support";
    public static final String MIN_NORMAL_COVERAGE_CFG = "min_normal_coverage";
    public static final String MIN_TUMOR_AF_CFG = "min_tumor_af";
    public static final String MAX_SHORT_STRAND_BIAS_CFG = "max_short_strand_bias";
    public static final String MIN_QUAL_BREAK_END_CFG = "min_qual_break_end";
    public static final String MIN_QUAL_BREAK_POINT_CFG = "min_qual_break_point";
    public static final String MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION = "min_qual_rescue_mobile_element_insertion";
    public static final String MAX_HOM_LENGTH_SHORT_INV_CFG = "max_hom_length_short_inv";
    public static final String MIN_LENGTH_CFG = "min_length";
    public static final String PON_DISTANCE = "pon_distance";
    private static final String FILTER_SGLS = "filter_sgls";

    public static FilterConstants from(final ConfigBuilder configBuilder)
    {
        RefGenomeVersion refGenVersion = RefGenomeVersion.from(configBuilder);

        return new FilterConstants(
                configBuilder.getInteger(HARD_MIN_TUMOR_QUAL_CFG),
                configBuilder.getInteger(HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_CFG),
                configBuilder.getDecimal(HARD_MAX_NORMAL_RELATIVE_SUPPORT_CFG),
                configBuilder.getDecimal(SOFT_MAX_NORMAL_RELATIVE_SUPPORT_CFG),
                configBuilder.getDecimal(MIN_NORMAL_COVERAGE_CFG),
                DEFAULT_MIN_TUMOR_AF_SGL,
                configBuilder.getDecimal(MIN_TUMOR_AF_CFG),
                configBuilder.getDecimal(MAX_SHORT_STRAND_BIAS_CFG),
                configBuilder.getInteger(MIN_QUAL_BREAK_END_CFG),
                configBuilder.getInteger(MIN_QUAL_BREAK_POINT_CFG),
                configBuilder.getInteger(MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION),
                configBuilder.getInteger(MAX_HOM_LENGTH_SHORT_INV_CFG),
                configBuilder.getInteger(MIN_LENGTH_CFG),
                configBuilder.getInteger(PON_DISTANCE),
                getPolyGRegions(refGenVersion), refGenVersion == V37 ? PMS2_V37 : PMS2_V38,
                configBuilder.hasFlag(FILTER_SGLS));
    }

    public FilterConstants(
            int minTumorQual, int hardMaxNormalAbsoluteSupport, double hardMaxNormalRelativeSupport, double softMaxNormalRelativeSupport,
            double minNormalCoverage, double minTumorAfBreakend, double minTumorAfBreakpoint, double maxShortStrandBias,
            int minQualBreakend, int minQualBreakpoint, int minQualRescueLine, int maxHomLengthShortInv,
            int minLength, int ponDistance, final List<ChrBaseRegion> polyGcRegions, final ChrBaseRegion lowQualRegion, boolean filterSGLs)
    {
        MinTumorQual = minTumorQual;
        HardMaxNormalAbsoluteSupport = hardMaxNormalAbsoluteSupport;
        HardMaxNormalRelativeSupport = hardMaxNormalRelativeSupport;
        SoftMaxNormalRelativeSupport = softMaxNormalRelativeSupport;
        MinNormalCoverage = minNormalCoverage;
        MinTumorAfBreakend = minTumorAfBreakend;
        MinTumorAfBreakpoint = minTumorAfBreakpoint;
        MaxShortStrandBias = maxShortStrandBias;
        MinQualBreakend = minQualBreakend;
        MinQualBreakpoint = minQualBreakpoint;
        MinQualRescueLine = minQualRescueLine;
        MaxHomLengthShortInv = maxHomLengthShortInv;
        MinLength = minLength;
        PonDistance = ponDistance;
        PolyGcRegions = polyGcRegions;
        LowQualRegion = lowQualRegion;
        FilterSGLs = filterSGLs;
    }

    public boolean matchesPolyGRegion(final String chromosome, int position)
    {
        return PolyGcRegions.stream().anyMatch(x -> x.containsPosition(chromosome, position));
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(HARD_MIN_TUMOR_QUAL_CFG, "Hard min tumor qual", DEFAULT_HARD_MIN_TUMOR_QUAL);

        configBuilder.addInteger(
                HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_CFG, "Hard max normal absolute support", DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT);

        configBuilder.addDecimal(
                HARD_MAX_NORMAL_RELATIVE_SUPPORT_CFG, "Hard max normal relative support", DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT);

        configBuilder.addDecimal(
                SOFT_MAX_NORMAL_RELATIVE_SUPPORT_CFG, "Soft max normal relative support", DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT);

        configBuilder.addDecimal(MIN_NORMAL_COVERAGE_CFG, "Min normal coverage", DEFAULT_MIN_NORMAL_COVERAGE);
        configBuilder.addDecimal(MIN_TUMOR_AF_CFG, "Min tumor allelic frequency for non-SGLs", DEFAULT_MIN_TUMOR_AF);
        configBuilder.addDecimal(MAX_SHORT_STRAND_BIAS_CFG, "Max short strand bias", DEFAULT_MAX_SHORT_STRAND_BIAS);
        configBuilder.addInteger(MIN_QUAL_BREAK_END_CFG, "Min qual break end", DEFAULT_MIN_QUAL_BREAK_END);
        configBuilder.addInteger(MIN_QUAL_BREAK_POINT_CFG, "Min qual break point", DEFAULT_MIN_QUAL_BREAK_POINT);

        configBuilder.addInteger(
                MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, "Min qual rescue mobile element insertions",
                DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION);

        configBuilder.addInteger(
                MAX_HOM_LENGTH_SHORT_INV_CFG, "Max homology length short inversion", DEFAULT_MAX_HOM_LENGTH_SHORT_INV);

        configBuilder.addInteger(MIN_LENGTH_CFG, "Min length", DEFAULT_MIN_LENGTH);
        configBuilder.addInteger(PON_DISTANCE, "PON permitted margin", DEFAULT_PON_DISTANCE);
        configBuilder.addFlag(FILTER_SGLS, "Filter SGLs from VCF, intended for tumor-only mode");
    }
}
