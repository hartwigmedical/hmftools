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
    String PLOIDY_PENALTY_STANDARD_DEVIATION = "ploidy_penalty_standard_deviation";
    String PLOIDY_PENALTY_MIN_STANDARD_DEVIATION = "ploidy_penalty_min_standard_deviation_per_ploidy";
    String PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER = "ploidy_penalty_major_allele_sub_one_multiplier";
    String PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL = "ploidy_penalty_major_allele_sub_one_additional";
    String PLOIDY_PENALTY_BASELINE_DEVIATION = "ploidy_penalty_baseline_deviation";

    double PLOIDY_PENALTY_FACTOR_DEFAULT = 0.3;
    double PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT = 0.05;
    double PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT = 1;
    double PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER_DEFAULT = 2;
    double PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL_DEFAULT = 1.5;
    double PLOIDY_PENALTY_BASELINE_DEVIATION_DEFAULT = 0.2;

    static void addOptions(@NotNull Options options) {
        options.addOption(PLOIDY_PENALTY_FACTOR, true, "Ploidy Penalty Factor");
        options.addOption(PLOIDY_PENALTY_STANDARD_DEVIATION, true, "Ploidy Penalty Standard Deviation");
        options.addOption(PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, true, "Ploidy Penalty Min Standard Deviation Per Ploidy");
        options.addOption(PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER, true, "PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER");
        options.addOption(PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL, true, "PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL");
        options.addOption(PLOIDY_PENALTY_BASELINE_DEVIATION, true, "PLOIDY_PENALTY_BASELINE_DEVIATION");
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
        final double ploidyPenaltyStandardDevation = defaultValue(cmd, PLOIDY_PENALTY_STANDARD_DEVIATION, PLOIDY_PENALTY_STANDARD_DEVIATION_DEFAULT);
        final double ploidyPenaltyMinStandardDevationPerPloidy = defaultValue(cmd, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION, PLOIDY_PENALTY_MIN_STANDARD_DEVIATION_DEFAULT);
        final double majorAlleleSubOnePenaltyMultiplier = defaultValue(cmd, PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER, PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_MULTIPLIER_DEFAULT);
        final double majorAlleleSubOneAdditionalPenalty = defaultValue(cmd, PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL, PLOIDY_PENALTY_MAJOR_ALLELE_SUB_ONE_ADDITIONAL_DEFAULT);
        final double baselineDeviation = defaultValue(cmd, PLOIDY_PENALTY_BASELINE_DEVIATION, PLOIDY_PENALTY_BASELINE_DEVIATION_DEFAULT);

        return ImmutableFitScoreConfig.builder()
                .ploidyPenaltyFactor(ploidyPenaltyFactor)
                .ploidyPenaltyStandardDeviation(ploidyPenaltyStandardDevation)
                .ploidyPenaltyMinStandardDeviationPerPloidy(ploidyPenaltyMinStandardDevationPerPloidy)
                .ploidyPenaltyMajorAlleleSubOneMultiplier(majorAlleleSubOnePenaltyMultiplier)
                .ploidyPenaltyMajorAlleleSubOneAdditional(majorAlleleSubOneAdditionalPenalty)
                .ploidyPenaltyBaselineDeviation(baselineDeviation)
                .build();
    }
}
