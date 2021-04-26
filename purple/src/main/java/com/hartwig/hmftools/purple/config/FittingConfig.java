package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

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

    private static final String MIN_PURITY = "min_purity";
    private static final String MAX_PURITY = "max_purity";
    private static final String PURITY_INCREMENT = "purity_increment";
    private static final String MIN_PLOIDY = "min_ploidy";
    private static final String MAX_PLOIDY = "max_ploidy";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";

    public static final double MIN_PURITY_DEFAULT = 0.08;
    public static final double MAX_PURITY_DEFAULT = 1.0;
    public static final double PURITY_INCREMENT_DEFAULT = 0.01;
    public static final double MIN_PLOIDY_DEFAULT = 1.0;
    public static final double MAX_PLOIDY_DEFAULT = 8;
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 150;

    // fitting scores
    private static final String PLOIDY_PENALTY_FACTOR = "ploidy_penalty_factor";
    private static final String PLOIDY_PENALTY_SUB_MIN_ADDITIONAL = "ploidy_penalty_sub_min_additional";
    private static final String PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER = "ploidy_penalty_sub_one_major_allele_multiplier";
    private static final String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    private static final String PLOIDY_PENALTY_MIN_STANDARD_DEVIATION = "ploidy_penalty_min_standard_deviation_per_ploidy";
    private static final String PLOIDY_PENALTY_MIN = "ploidy_penalty_min";

    double PLOIDY_PENALTY_FACTOR_DEFAULT = 0.4;
    double PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.05;
    double PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT = 1.5;
    double PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT = 1;
    double PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT = 1.5;
    double PLOIDY_PENALTY_MIN_DEFAULT = 0.1;

    public FittingConfig(@NotNull final CommandLine cmd)
    {
        MinPurity = getConfigValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
        MaxPurity = getConfigValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
        PurityIncrement = getConfigValue(cmd, PURITY_INCREMENT, PURITY_INCREMENT_DEFAULT);
        MinPloidy = getConfigValue(cmd, MIN_PLOIDY, MIN_PLOIDY_DEFAULT);
        MaxPloidy = getConfigValue(cmd, MAX_PLOIDY, MAX_PLOIDY_DEFAULT);

        MinDiploidTumorRatioCount = getConfigValue(cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT, MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT);
        MinDiploidTumorRatioCountAtCentromere = getConfigValue(
                cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE, MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT);

        PloidyPenaltyFactor = getConfigValue(cmd, PLOIDY_PENALTY_FACTOR, PLOIDY_PENALTY_FACTOR_DEFAULT);

        PloidyPenaltyStandardDeviation =
                getConfigValue(cmd, PLOIDY_PENALTY_STANDARD_DEVIATION, PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);

        PloidyPenaltyMinStandardDeviationPerPloidy =
                getConfigValue(cmd, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT);

        PloidyPenaltyMajorAlleleSubOneMultiplier =
                getConfigValue(cmd, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT);

        PloidyPenaltyMajorAlleleSubOneAdditional =
                getConfigValue(cmd, PLOIDY_PENALTY_SUB_MIN_ADDITIONAL, PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT);

        PloidyPenaltyBaselineDeviation = getConfigValue(cmd, PLOIDY_PENALTY_MIN, PLOIDY_PENALTY_MIN_DEFAULT);
    }

    static void addOptions(@NotNull Options options)
    {
        options.addOption(MIN_PURITY, true, "Minimum purity (default " + MIN_PURITY_DEFAULT + ")");
        options.addOption(MAX_PURITY, true, "Maximum purity (default " + MAX_PURITY_DEFAULT + ")");
        options.addOption(PURITY_INCREMENT, true, "Purity increment (default " + PURITY_INCREMENT_DEFAULT + ")");

        options.addOption(MIN_PLOIDY, true, "Minimum ploidy (default " + MIN_PLOIDY_DEFAULT + ")");
        options.addOption(MAX_PLOIDY, true, "Maximum ploidy (default " + MAX_PLOIDY_DEFAULT + ")");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere.");

        options.addOption(PLOIDY_PENALTY_FACTOR, true, "Penalty factor to apply to the number of copy number events");
        options.addOption(PLOIDY_PENALTY_STANDARD_DEVIATION,
                true,
                "Standard deviation of normal distribution modelling ploidy deviation from whole number");
        options.addOption(PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, true, "Minimum ploidy penalty standard deviation to be applied");
        options.addOption(PLOIDY_PENALTY_SUB_MIN_ADDITIONAL, true, "Additional penalty to apply to major allele < 1 or minor allele < 0");
        options.addOption(PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER, true, "Penalty multiplier applied to major allele < 1 ");
        options.addOption(PLOIDY_PENALTY_MIN, true, "Minimum ploidy penalty");
    }
}
