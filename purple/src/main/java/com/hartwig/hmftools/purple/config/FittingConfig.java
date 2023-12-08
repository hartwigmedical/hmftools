package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.config.PurpleConstants.DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE;
import static com.hartwig.hmftools.purple.config.PurpleConstants.DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MAX_PLOIDY_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MAX_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MIN_PLOIDY_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.PURITY_INCREMENT_DEFAULT;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FittingConfig
{
    public final double MinPurity;
    public final double MaxPurity;
    public final double PurityIncrement;
    public final double MinPloidy;
    public final double MaxPloidy;

    public final int MinDiploidTumorRatioCount;
    public final int MinDiploidTumorRatioCountAtCentromere;

    public final double PloidyPenaltyFactor;
    public final double PloidyPenaltyStandardDeviation;
    public final double PloidyPenaltyMinStandardDeviationPerPloidy;
    public final double PloidyPenaltyMajorAlleleSubOneMultiplier;
    public final double PloidyPenaltyMajorAlleleSubOneAdditional;
    public final double PloidyPenaltyBaselineDeviation;

    public final int RecoveryMinMateQualScore;
    public final int RecoveryMinSglQualScore;

    public final double DeviationPenaltyGcAdjust;
    public final double GcRatioExponent;

    private static final String MIN_PURITY = "min_purity";
    private static final String MAX_PURITY = "max_purity";
    private static final String PURITY_INCREMENT = "purity_increment";
    private static final String MIN_PLOIDY = "min_ploidy";
    private static final String MAX_PLOIDY = "max_ploidy";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";

    private static final String DEVIATION_PENALTY_GC_ADJUST = "deviation_penalty_gc_adjust";
    private static final String GC_RATIO_EXPONENT = "gc_ratio_exponent";

    // fitting scores
    private static final String PLOIDY_PENALTY_FACTOR = "ploidy_penalty_factor";
    private static final String PLOIDY_PENALTY_SUB_MIN_ADDITIONAL = "ploidy_penalty_sub_min_additional";
    private static final String PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER = "ploidy_penalty_sub_one_major_allele_multiplier";
    private static final String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    private static final String PLOIDY_PENALTY_MIN_STANDARD_DEVIATION = "ploidy_penalty_min_standard_deviation_per_ploidy";
    private static final String PLOIDY_PENALTY_MIN = "ploidy_penalty_min";

    private static final double PLOIDY_PENALTY_FACTOR_DEFAULT = 0.4;
    private static final double PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.05;
    private static final double PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT = 1.5;
    private static final double PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT = 1;
    private static final double PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT = 1.5;
    private static final double PLOIDY_PENALTY_MIN_DEFAULT = 0.1;

    // SV recovery
    public static final String CFG_MIN_MATE_QUAL_SCORE = "recovery_mate_min_qual";
    public static final String CFG_MIN_SGL_QUAL_SCORE = "recovery_sgl_min_qual";

    public FittingConfig(final ConfigBuilder configBuilder)
    {
        MinPurity = configBuilder.getDecimal(MIN_PURITY);
        MaxPurity = configBuilder.getDecimal(MAX_PURITY);
        PurityIncrement = configBuilder.getDecimal(PURITY_INCREMENT);
        MinPloidy = configBuilder.getDecimal(MIN_PLOIDY);
        MaxPloidy = configBuilder.getDecimal(MAX_PLOIDY);

        MinDiploidTumorRatioCount = configBuilder.getInteger(MIN_DIPLOID_TUMOR_RATIO_COUNT);

        MinDiploidTumorRatioCountAtCentromere = configBuilder.getInteger(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE);

        PloidyPenaltyFactor = configBuilder.getDecimal(PLOIDY_PENALTY_FACTOR);

        PloidyPenaltyStandardDeviation = configBuilder.getDecimal(PLOIDY_PENALTY_STANDARD_DEVIATION);
        PloidyPenaltyMinStandardDeviationPerPloidy = configBuilder.getDecimal(PLOIDY_PENALTY_MIN_STANDARD_DEVIATION);
        PloidyPenaltyMajorAlleleSubOneMultiplier = configBuilder.getDecimal(PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER);
        PloidyPenaltyMajorAlleleSubOneAdditional = configBuilder.getDecimal(PLOIDY_PENALTY_SUB_MIN_ADDITIONAL);

        PloidyPenaltyBaselineDeviation = configBuilder.getDecimal(PLOIDY_PENALTY_MIN);

        RecoveryMinMateQualScore = configBuilder.getInteger(CFG_MIN_MATE_QUAL_SCORE);
        RecoveryMinSglQualScore = configBuilder.getInteger(CFG_MIN_SGL_QUAL_SCORE);

        DeviationPenaltyGcAdjust = configBuilder.getDecimal(DEVIATION_PENALTY_GC_ADJUST);
        GcRatioExponent = configBuilder.getDecimal(GC_RATIO_EXPONENT);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(MIN_PURITY, "Minimum purity", MIN_PURITY_DEFAULT);
        configBuilder.addDecimal(MAX_PURITY, "Maximum purity", MAX_PURITY_DEFAULT);
        configBuilder.addDecimal(PURITY_INCREMENT, "Purity increment", PURITY_INCREMENT_DEFAULT);
        configBuilder.addDecimal(MIN_PLOIDY, "Minimum ploidy", MIN_PLOIDY_DEFAULT);
        configBuilder.addDecimal(MAX_PLOIDY, "Maximum ploidy", MAX_PLOIDY_DEFAULT);

        configBuilder.addInteger(
                MIN_DIPLOID_TUMOR_RATIO_COUNT,
                "Minimum ratio count while smoothing before diploid regions become suspect",
                MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT);

        configBuilder.addInteger(
                MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere",
                MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_FACTOR, "Penalty factor to apply to the number of copy number events", PLOIDY_PENALTY_FACTOR_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_STANDARD_DEVIATION,
                "Standard deviation of normal distribution modelling ploidy deviation from whole number",
                PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_MIN_STANDARD_DEVIATION,
                "Minimum ploidy penalty standard deviation to be applied", PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_SUB_MIN_ADDITIONAL,
                "Additional penalty to apply to major allele < 1 or minor allele < 0",
                PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER,
                "Penalty multiplier applied to major allele < 1", PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT);

        configBuilder.addDecimal(PLOIDY_PENALTY_MIN, "Minimum ploidy penalty", PLOIDY_PENALTY_MIN_DEFAULT);

        configBuilder.addInteger(
                CFG_MIN_MATE_QUAL_SCORE, "SV recovery non-SGL min qual score", DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE);

        configBuilder.addInteger(
                CFG_MIN_SGL_QUAL_SCORE, "SV recovery SGL min qual score", DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE);

        configBuilder.addDecimal(DEVIATION_PENALTY_GC_ADJUST, "Adjust deviation penalty by tumor GC Ratio", 0);
        configBuilder.addDecimal(GC_RATIO_EXPONENT, "Adjust GC Ratio by exponent in penalty calc", 0);
    }
}
