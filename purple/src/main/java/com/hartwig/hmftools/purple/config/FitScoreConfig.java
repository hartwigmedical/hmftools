package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FitScoreConfig {

    String PLOIDY_PENALTY_FACTOR = "ploidy_penalty_factor";
    String PLOIDY_PENALTY_SUB_MIN_ADDITIONAL = "ploidy_penalty_sub_min_additional";
    String PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER = "ploidy_penalty_sub_one_major_allele_multiplier";
    String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    String PLOIDY_PENALTY_MIN_STANDARD_DEVIATION = "ploidy_penalty_min_standard_deviation_per_ploidy";
    String PLOIDY_PENALTY_MIN = "ploidy_penalty_min";

    double PLOIDY_PENALTY_FACTOR_DEFAULT = 0.4;
    double PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.05;
    double PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT = 1.5;
    double PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT = 1;
    double PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT = 1.5;
    double PLOIDY_PENALTY_MIN_DEFAULT = 0.1;

    static void addOptions(@NotNull Options options) {
        options.addOption(PLOIDY_PENALTY_FACTOR, true, "Penalty factor to apply to the number of copy number events");
        options.addOption(PLOIDY_PENALTY_STANDARD_DEVIATION,
                true,
                "Standard deviation of normal distribution modelling ploidy deviation from whole number");
        options.addOption(PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, true, "Minimum ploidy penalty standard deviation to be applied");
        options.addOption(PLOIDY_PENALTY_SUB_MIN_ADDITIONAL, true, "Additional penalty to apply to major allele < 1 or minor allele < 0");
        options.addOption(PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER, true, "Penalty multiplier applied to major allele < 1 ");
        options.addOption(PLOIDY_PENALTY_MIN, true, "Minimum ploidy penalty");
    }

    double ploidyPenaltyFactor();

    double ploidyPenaltyStandardDeviation();

    double ploidyPenaltyMinStandardDeviationPerPloidy();

    double ploidyPenaltyMajorAlleleSubOneMultiplier();

    double ploidyPenaltyMajorAlleleSubOneAdditional();

    double ploidyPenaltyBaselineDeviation();

    @NotNull
    static FitScoreConfig createConfig(@NotNull final CommandLine cmd) {
        final double ploidyPenaltyFactor = defaultValue(cmd, PLOIDY_PENALTY_FACTOR, PLOIDY_PENALTY_FACTOR_DEFAULT);
        final double ploidyPenaltyStandardDeviation =
                defaultValue(cmd, PLOIDY_PENALTY_STANDARD_DEVIATION, PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);
        final double ploidyPenaltyMinStandardDeviationPerPloidy =
                defaultValue(cmd, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT);
        final double majorAlleleSubOnePenaltyMultiplier =
                defaultValue(cmd, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER, PLOIDY_PENALTY_SUB_ONE_MAJOR_ALLELE_MULTIPLIER_DEFAULT);
        final double majorAlleleSubOneAdditionalPenalty =
                defaultValue(cmd, PLOIDY_PENALTY_SUB_MIN_ADDITIONAL, PLOIDY_PENALTY_SUB_MIN_ADDITIONAL_DEFAULT);
        final double baselineDeviation = defaultValue(cmd, PLOIDY_PENALTY_MIN, PLOIDY_PENALTY_MIN_DEFAULT);

        return ImmutableFitScoreConfig.builder()
                .ploidyPenaltyFactor(ploidyPenaltyFactor)
                .ploidyPenaltyStandardDeviation(ploidyPenaltyStandardDeviation)
                .ploidyPenaltyMinStandardDeviationPerPloidy(ploidyPenaltyMinStandardDeviationPerPloidy)
                .ploidyPenaltyMajorAlleleSubOneMultiplier(majorAlleleSubOnePenaltyMultiplier)
                .ploidyPenaltyMajorAlleleSubOneAdditional(majorAlleleSubOneAdditionalPenalty)
                .ploidyPenaltyBaselineDeviation(baselineDeviation)
                .build();
    }
}
