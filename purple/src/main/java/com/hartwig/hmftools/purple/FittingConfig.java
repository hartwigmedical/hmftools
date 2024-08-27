package com.hartwig.hmftools.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigBuilder.getConfigDecimal;
import static com.hartwig.hmftools.common.utils.config.ConfigBuilder.getConfigInteger;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.purple.PurpleConstants.DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE;
import static com.hartwig.hmftools.purple.PurpleConstants.DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE;
import static com.hartwig.hmftools.purple.PurpleConstants.MAX_PLOIDY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.MAX_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_PLOIDY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_FACTOR_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_MIN_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.PURITY_INCREMENT_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_GC_RATIO_EXPONENT_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_PLOIDY_PENALTY_MIN_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT;

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
    public final double PloidyPenaltyMinDeviation;

    public final int RecoveryMinMateQualScore;
    public final int RecoveryMinSglQualScore;

    public final double DeviationPenaltyGcMinAdjust;
    public final double GcRatioExponent;

    public static final String MIN_PURITY = "min_purity";
    public static final String MAX_PURITY = "max_purity";
    public static final String PURITY_INCREMENT = "purity_increment";
    public static final String MIN_PLOIDY = "min_ploidy";
    public static final String MAX_PLOIDY = "max_ploidy";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";

    private static final String DEVIATION_PENALTY_GC_MIN_ADJUST = "deviation_penalty_gc_min_adjust";
    private static final String GC_RATIO_EXPONENT = "gc_ratio_exponent";

    // fitting scores
    private static final String PLOIDY_PENALTY_FACTOR = "ploidy_penalty_factor";
    private static final String PLOIDY_PENALTY_SUB_MIN_ADDITIONAL = "ploidy_penalty_sub_min_additional";
    private static final String PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER = "ploidy_penalty_sub_one_major_allele_multiplier";
    private static final String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    private static final String PLOIDY_PENALTY_MIN_STANDARD_DEVIATION = "ploidy_penalty_min_standard_deviation_per_ploidy";
    private static final String PLOIDY_PENALTY_MIN = "ploidy_penalty_min";

    // SV recovery
    public static final String CFG_MIN_MATE_QUAL_SCORE = "recovery_mate_min_qual";
    public static final String CFG_MIN_SGL_QUAL_SCORE = "recovery_sgl_min_qual";

    public FittingConfig(final ConfigBuilder configBuilder, boolean targetedMode)
    {
        MinPurity = configBuilder.getDecimal(MIN_PURITY);
        MaxPurity = configBuilder.getDecimal(MAX_PURITY);
        PurityIncrement = configBuilder.getDecimal(PURITY_INCREMENT);
        MinPloidy = configBuilder.getDecimal(MIN_PLOIDY);
        MaxPloidy = configBuilder.getDecimal(MAX_PLOIDY);

        MinDiploidTumorRatioCount = getConfigInteger(
                configBuilder, MIN_DIPLOID_TUMOR_RATIO_COUNT,
                targetedMode ? TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT : MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT);

        MinDiploidTumorRatioCountAtCentromere = getConfigInteger(
                configBuilder, MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                targetedMode ? TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT : MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT);

        PloidyPenaltyFactor = configBuilder.getDecimal(PLOIDY_PENALTY_FACTOR);

        PloidyPenaltyStandardDeviation = getConfigDecimal(
                configBuilder, PLOIDY_PENALTY_STANDARD_DEVIATION,
                targetedMode ? TARGETED_PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT : PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);

        PloidyPenaltyMinStandardDeviationPerPloidy = configBuilder.getDecimal(PLOIDY_PENALTY_MIN_STANDARD_DEVIATION);

        PloidyPenaltyMajorAlleleSubOneMultiplier = getConfigDecimal(
                configBuilder, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER,
                targetedMode ? TARGETED_PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT : PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT);

        PloidyPenaltyMajorAlleleSubOneAdditional = configBuilder.getDecimal(PLOIDY_PENALTY_SUB_MIN_ADDITIONAL);

        PloidyPenaltyMinDeviation = getConfigDecimal(
                configBuilder, PLOIDY_PENALTY_MIN, targetedMode ? TARGETED_PLOIDY_PENALTY_MIN_DEFAULT : PLOIDY_PENALTY_MIN_DEFAULT);

        RecoveryMinMateQualScore = configBuilder.getInteger(CFG_MIN_MATE_QUAL_SCORE);
        RecoveryMinSglQualScore = configBuilder.getInteger(CFG_MIN_SGL_QUAL_SCORE);

        DeviationPenaltyGcMinAdjust = getConfigDecimal(
                configBuilder, DEVIATION_PENALTY_GC_MIN_ADJUST, targetedMode ? TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT : 0);

        GcRatioExponent = getConfigDecimal(configBuilder, GC_RATIO_EXPONENT, targetedMode ? TARGETED_GC_RATIO_EXPONENT_DEFAULT : 0);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(MIN_PURITY, "Minimum purity", MIN_PURITY_DEFAULT);
        configBuilder.addDecimal(MAX_PURITY, "Maximum purity", MAX_PURITY_DEFAULT);
        configBuilder.addDecimal(PURITY_INCREMENT, "Purity increment", PURITY_INCREMENT_DEFAULT);
        configBuilder.addDecimal(MIN_PLOIDY, "Minimum ploidy", MIN_PLOIDY_DEFAULT);
        configBuilder.addDecimal(MAX_PLOIDY, "Maximum ploidy", MAX_PLOIDY_DEFAULT);

        addTargetedInteger(
                configBuilder, MIN_DIPLOID_TUMOR_RATIO_COUNT,
                "Minimum ratio count while smoothing before diploid regions become suspect",
                MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT, TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT);

        addTargetedInteger(
                configBuilder, MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere",
                MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT, TARGETED_MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_FACTOR, "Penalty factor to apply to the number of copy number events", PLOIDY_PENALTY_FACTOR_DEFAULT);

        addTargetedDecimal(
                configBuilder, PLOIDY_PENALTY_STANDARD_DEVIATION,
                "Standard deviation of normal distribution modelling ploidy deviation from whole number",
                PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT, TARGETED_PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_MIN_STANDARD_DEVIATION,
                "Minimum ploidy penalty standard deviation to be applied", PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT);

        configBuilder.addDecimal(
                PLOIDY_PENALTY_SUB_MIN_ADDITIONAL,
                "Additional penalty to apply to major allele < 1 or minor allele < 0",
                PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT);

        addTargetedDecimal(
                configBuilder, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER, "Penalty multiplier applied to major allele < 1",
                PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT, TARGETED_PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT);

        addTargetedDecimal(
                configBuilder, PLOIDY_PENALTY_MIN, "Minimum ploidy penalty",
                PLOIDY_PENALTY_MIN_DEFAULT, TARGETED_PLOIDY_PENALTY_MIN_DEFAULT);

        configBuilder.addInteger(
                CFG_MIN_MATE_QUAL_SCORE, "SV recovery non-SGL min qual score", DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE);

        configBuilder.addInteger(
                CFG_MIN_SGL_QUAL_SCORE, "SV recovery SGL min qual score", DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE);

        addTargetedDecimal(
                configBuilder, DEVIATION_PENALTY_GC_MIN_ADJUST, "Adjust deviation penalty by tumor GC Ratio",
                0, TARGETED_DEVIATION_PENALTY_GC_MIN_ADJUST_DEFAULT);

        addTargetedDecimal(
                configBuilder, GC_RATIO_EXPONENT, "Adjust GC Ratio by exponent in penalty calc",
                0, TARGETED_GC_RATIO_EXPONENT_DEFAULT);
    }

    private static void addTargetedDecimal(
            final ConfigBuilder configBuilder, final String name, final String desc, double defaultValue, double targetedDefaultValue)
    {
        configBuilder.addConfigItem(
                DECIMAL, name, false,
                format("%s, default=%.3g targeted default=%.3g", desc, defaultValue, targetedDefaultValue), String.valueOf(defaultValue));
    }

    private static void addTargetedInteger(
            final ConfigBuilder configBuilder, final String name, final String desc, int defaultValue, int targetedDefaultValue)
    {
        configBuilder.addConfigItem(
                INTEGER, name, false,
                format("%s, default=%d targeted default=%d", desc, defaultValue, targetedDefaultValue), String.valueOf(defaultValue));
    }
}
