package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;

public class FilterConstants
{
    public final int MinTumorQual;
    public final int HardMaxNormalAbsoluteSupport;
    public final double HardMaxNormalRelativeSupport;
    public final double SoftMaxNormalRelativeSupport;
    public final double MinNormalCoverage;
    public final double MinTumorAf;
    public final double MaxShortStrandBias;
    public final int MinQualBreakend;
    public final int MinQualBreakpoint;
    public final int MinQualRescueLine;
    public final int MaxHomLengthShortInv;
    public final int MaxInexactHomLengthShortDel;
    public final int MinLength;
    public final int PonDistance;
    public final ChrBaseRegion PolyGcRegion;

    // filter which only apply when reference is present:
    // minNormalCoverage, minRelativeCoverage, maxNormalSupport, shortSRNormalSupport, discordantPairSupport

    // default filter values
    public static final int SHORT_RESCUE_LENGTH = 1000;
    public static final int SHORT_CALLING_SIZE = 1000;
    public static final int HOM_INV_LENGTH = 40;

    public static final String POLY_G_INSERT = "GGGGGGGGGGGGGGGG";
    public static final String POLY_C_INSERT = "CCCCCCCCCCCCCCCC";
    public static final String POLY_A_HOMOLOGY = "AAAAAAA";
    public static final String POLY_T_HOMOLOGY = "TTTTTTT";

    public static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 100;
    public static final int DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT = 3;
    public static final double DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT = 0.08;
    public static final double DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT = 0.03;
    public static final double DEFAULT_MIN_NORMAL_COVERAGE = 8;
    public static final double DEFAULT_MIN_TUMOR_AF = 0.005;
    public static final double DEFAULT_MAX_SHORT_STRAND_BIAS = 0.95;
    public static final int DEFAULT_MIN_QUAL_BREAK_END = 1000;
    public static final int DEFAULT_MIN_QUAL_BREAK_POINT = 400;
    public static final int DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION = 1000;
    public static final int DEFAULT_MAX_HOM_LENGTH_SHORT_INV = 6;
    public static final int DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL = 5;
    public static final int DEFAULT_MIN_LENGTH = 32;
    public static final int DEFAULT_PON_DISTANCE = 3;

    public static final String POLY_A = "AAAAAAA";
    public static final String POLY_T = "TTTTTTT";

    public static final ChrBaseRegion LINC_00486_V37 = new ChrBaseRegion("2", 33141260, 33141700);
    public static final ChrBaseRegion LINC_00486_V38 = new ChrBaseRegion("chr2", 32916190, 32916630);

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
    public static final String MAX_INEXACT_HOM_LENGTH_SHORT_DEL_CFG = "max_inexact_hom_length_short_del";
    public static final String MIN_LENGTH_CFG = "min_length";
    public static final String PON_DISTANCE = "pon_distance";

    public static FilterConstants from(final CommandLine cmd)
    {
        RefGenomeVersion refGenVersion = RefGenomeVersion.valueOf(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        return new FilterConstants(
                Integer.parseInt(cmd.getOptionValue(HARD_MIN_TUMOR_QUAL_CFG, String.valueOf(DEFAULT_HARD_MIN_TUMOR_QUAL))),
                Integer.parseInt(cmd.getOptionValue(
                        HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_CFG, String.valueOf(DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT))),
                Double.parseDouble(cmd.getOptionValue(
                        HARD_MAX_NORMAL_RELATIVE_SUPPORT_CFG, String.valueOf(DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT))),
                Double.parseDouble(cmd.getOptionValue(
                        SOFT_MAX_NORMAL_RELATIVE_SUPPORT_CFG, String.valueOf(DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT))),
                Double.parseDouble(cmd.getOptionValue(MIN_NORMAL_COVERAGE_CFG, String.valueOf(DEFAULT_MIN_NORMAL_COVERAGE))),
                Double.parseDouble(cmd.getOptionValue(MIN_TUMOR_AF_CFG, String.valueOf(DEFAULT_MIN_TUMOR_AF))),
                Double.parseDouble(cmd.getOptionValue(MAX_SHORT_STRAND_BIAS_CFG, String.valueOf(DEFAULT_MAX_SHORT_STRAND_BIAS))),
                Integer.parseInt(cmd.getOptionValue(MIN_QUAL_BREAK_END_CFG, String.valueOf(DEFAULT_MIN_QUAL_BREAK_END))),
                Integer.parseInt(cmd.getOptionValue(MIN_QUAL_BREAK_POINT_CFG, String.valueOf(DEFAULT_MIN_QUAL_BREAK_POINT))),
                Integer.parseInt(cmd.getOptionValue(
                        MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, String.valueOf(DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION))),
                Integer.parseInt(cmd.getOptionValue(MAX_HOM_LENGTH_SHORT_INV_CFG, String.valueOf(DEFAULT_MAX_HOM_LENGTH_SHORT_INV))),
                Integer.parseInt(cmd.getOptionValue(
                        MAX_INEXACT_HOM_LENGTH_SHORT_DEL_CFG, String.valueOf(DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL))),
                Integer.parseInt(cmd.getOptionValue(MIN_LENGTH_CFG, String.valueOf(DEFAULT_MIN_LENGTH))),
                Integer.parseInt(cmd.getOptionValue(PON_DISTANCE, String.valueOf(DEFAULT_PON_DISTANCE))),
                refGenVersion == V37 ? LINC_00486_V37 : LINC_00486_V38);
    }

    public FilterConstants(
            int minTumorQual, int hardMaxNormalAbsoluteSupport, double hardMaxNormalRelativeSupport, double softMaxNormalRelativeSupport,
            double minNormalCoverage, double minTumorAf, double maxShortStrandBias, int minQualBreakend, int minQualBreakpoint,
            int minQualRescueLine, int maxHomLengthShortInv, int maxInexactHomLengthShortDel, int minLength, int ponDistance,
            final ChrBaseRegion polyGcRegion)
    {
        MinTumorQual = minTumorQual;
        HardMaxNormalAbsoluteSupport = hardMaxNormalAbsoluteSupport;
        HardMaxNormalRelativeSupport = hardMaxNormalRelativeSupport;
        SoftMaxNormalRelativeSupport = softMaxNormalRelativeSupport;
        MinNormalCoverage = minNormalCoverage;
        MinTumorAf = minTumorAf;
        MaxShortStrandBias = maxShortStrandBias;
        MinQualBreakend = minQualBreakend;
        MinQualBreakpoint = minQualBreakpoint;
        MinQualRescueLine = minQualRescueLine;
        MaxHomLengthShortInv = maxHomLengthShortInv;
        MaxInexactHomLengthShortDel = maxInexactHomLengthShortDel;
        MinLength = minLength;
        PonDistance = ponDistance;
        PolyGcRegion = polyGcRegion;
    }
}
