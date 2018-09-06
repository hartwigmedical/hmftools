package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FittingConfig {

    String MIN_PURITY = "min_purity";
    String MAX_PURITY = "max_purity";
    String PURITY_INCREMENT = "purity_increment";
    String MIN_NORM_FACTOR = "min_norm_factor";
    String MAX_NORM_FACTOR = "max_norm_factor";
    String NORM_FACTOR_INCREMENTS = "norm_factor_increment";

    double MIN_PURITY_DEFAULT = 0.08;
    double MAX_PURITY_DEFAULT = 1.0;
    double PURITY_INCREMENT_DEFAULT = 0.01;
    double MIN_NORM_FACTOR_DEFAULT = 0.33;
    double MAX_NORM_FACTOR_DEFAULT = 2.0;
    double NORM_FACTOR_INCREMENTS_DEFAULT = 0.01;

    static void addOptions(@NotNull Options options) {
        options.addOption(MIN_PURITY, true, "Minimum purity (default " + MIN_PURITY_DEFAULT + ")");
        options.addOption(MAX_PURITY, true, "Maximum purity (default " + MAX_PURITY_DEFAULT + ")");
        options.addOption(PURITY_INCREMENT, true, "Purity increment (default " + PURITY_INCREMENT_DEFAULT + ")");

        options.addOption(MIN_NORM_FACTOR, true, "Minimum norm factor (default " + MIN_NORM_FACTOR_DEFAULT + ")");
        options.addOption(MAX_NORM_FACTOR, true, "Maximum norm factor (default " + MAX_NORM_FACTOR_DEFAULT + ")");
        options.addOption(NORM_FACTOR_INCREMENTS, true, "Norm factor increments (default  " + NORM_FACTOR_INCREMENTS_DEFAULT + ")");
    }

    double minPurity();

    double maxPurity();

    double purityIncrement();

    double minNormFactor();

    double maxNormFactor();

    double normFactorIncrement();

    default int maxPloidy() {
        return 20;
    }

    @NotNull
    static FittingConfig createConfig(@NotNull final CommandLine cmd) {
        final double minPurity = defaultValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
        final double maxPurity = defaultValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
        final double purityIncrement = defaultValue(cmd, PURITY_INCREMENT, PURITY_INCREMENT_DEFAULT);
        final double minNormFactor = defaultValue(cmd, MIN_NORM_FACTOR, MIN_NORM_FACTOR_DEFAULT);
        final double maxNormFactor = defaultValue(cmd, MAX_NORM_FACTOR, MAX_NORM_FACTOR_DEFAULT);
        final double normFactorIncrement = defaultValue(cmd, NORM_FACTOR_INCREMENTS, NORM_FACTOR_INCREMENTS_DEFAULT);

        return ImmutableFittingConfig.builder()
                .minPurity(minPurity)
                .maxPurity(maxPurity)
                .purityIncrement(purityIncrement)
                .minNormFactor(minNormFactor)
                .maxNormFactor(maxNormFactor)
                .normFactorIncrement(normFactorIncrement)
                .build();

    }

}
