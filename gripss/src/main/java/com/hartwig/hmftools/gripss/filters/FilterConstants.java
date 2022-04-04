package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

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
    public static final int HOM_INV_LENGTH = 40;

    public static final int LINE_POLY_AT_TEST_LEN = 18;
    public static final int LINE_POLY_AT_REQ = 16;
    public static final int SGL_INS_SEQ_MIN_LENGTH = 16;
    public static final double SGL_MIN_STRAND_BIAS = 0.05;
    public static final double SGL_MAX_STRAND_BIAS = 0.95;

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
    public static final double DEFAULT_MIN_TUMOR_AF_SGL = 0.015;
    public static final double DEFAULT_MAX_SHORT_STRAND_BIAS = 0.95;
    public static final int DEFAULT_MIN_QUAL_BREAK_END = 500;
    public static final int DEFAULT_MIN_QUAL_BREAK_POINT = 400;
    public static final int DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION = 500;
    public static final int DEFAULT_MAX_HOM_LENGTH_SHORT_INV = 6;
    public static final int DEFAULT_MIN_LENGTH = 32;
    public static final int DEFAULT_PON_DISTANCE = 3;

    public static final String POLY_A = "AAAAAAA";
    public static final String POLY_T = "TTTTTTT";

    public static final List<ChrBaseRegion> POLY_G_REGIONS_V37 = Lists.newArrayList();
    public static final List<ChrBaseRegion> POLY_G_REGIONS_V38 = Lists.newArrayList();

    static
    {
        POLY_G_REGIONS_V37.add(new ChrBaseRegion("2", 33141260, 33141700));
        POLY_G_REGIONS_V37.add(new ChrBaseRegion("4", 41218427, 41218467));
        POLY_G_REGIONS_V37.add(new ChrBaseRegion("17", 42646418, 42646458));

        POLY_G_REGIONS_V38.add(new ChrBaseRegion("chr2", 32916190, 32916630));
        POLY_G_REGIONS_V38.add(new ChrBaseRegion("chr4", 41216410, 41216450));
        POLY_G_REGIONS_V38.add(new ChrBaseRegion("chr17", 44569050, 44569090));
    }

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
                DEFAULT_MIN_TUMOR_AF_SGL, Double.parseDouble(cmd.getOptionValue(MIN_TUMOR_AF_CFG, String.valueOf(DEFAULT_MIN_TUMOR_AF))),
                Double.parseDouble(cmd.getOptionValue(MAX_SHORT_STRAND_BIAS_CFG, String.valueOf(DEFAULT_MAX_SHORT_STRAND_BIAS))),
                Integer.parseInt(cmd.getOptionValue(MIN_QUAL_BREAK_END_CFG, String.valueOf(DEFAULT_MIN_QUAL_BREAK_END))),
                Integer.parseInt(cmd.getOptionValue(MIN_QUAL_BREAK_POINT_CFG, String.valueOf(DEFAULT_MIN_QUAL_BREAK_POINT))),
                Integer.parseInt(cmd.getOptionValue(
                        MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, String.valueOf(DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION))),
                Integer.parseInt(cmd.getOptionValue(MAX_HOM_LENGTH_SHORT_INV_CFG, String.valueOf(DEFAULT_MAX_HOM_LENGTH_SHORT_INV))),
                Integer.parseInt(cmd.getOptionValue(MIN_LENGTH_CFG, String.valueOf(DEFAULT_MIN_LENGTH))),
                Integer.parseInt(cmd.getOptionValue(PON_DISTANCE, String.valueOf(DEFAULT_PON_DISTANCE))),
                refGenVersion == V37 ? POLY_G_REGIONS_V37 : POLY_G_REGIONS_V38, refGenVersion == V37 ? PMS2_V37 : PMS2_V38,
                cmd.hasOption(FILTER_SGLS));
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

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(HARD_MIN_TUMOR_QUAL_CFG, true, "Hard min tumor qual");
        options.addOption(HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_CFG, true, "Hard max normal absolute support");
        options.addOption(HARD_MAX_NORMAL_RELATIVE_SUPPORT_CFG, true, "Hard max normal relative support");
        options.addOption(SOFT_MAX_NORMAL_RELATIVE_SUPPORT_CFG, true, "Soft max normal relative support");
        options.addOption(MIN_NORMAL_COVERAGE_CFG, true, "Min normal coverage");
        options.addOption(MIN_TUMOR_AF_CFG, true, "Min tumor allelic frequency for non-SGLs");
        options.addOption(MAX_SHORT_STRAND_BIAS_CFG, true, "Max short strand bias");
        options.addOption(MIN_QUAL_BREAK_END_CFG, true, "Min qual break end");
        options.addOption(MIN_QUAL_BREAK_POINT_CFG, true, "Min qual break point");
        options.addOption(MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, true, "Min qual rescue mobile element insertions");
        options.addOption(MAX_HOM_LENGTH_SHORT_INV_CFG, true, "Max homology length short inversion");
        options.addOption(MIN_LENGTH_CFG, true, "Min length");
        options.addOption(PON_DISTANCE, true, "PON permitted margin");
        options.addOption(FILTER_SGLS, false, "Filter SGLs from VCF, intended for tumor-only mode");
    }

}
