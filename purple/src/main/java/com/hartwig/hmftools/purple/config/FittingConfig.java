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
    String MIN_PLOIDY = "min_ploidy";
    String MAX_PLOIDY = "max_ploidy";

    double MIN_PURITY_DEFAULT = 0.08;
    double MAX_PURITY_DEFAULT = 1.0;
    double PURITY_INCREMENT_DEFAULT = 0.01;
    double MIN_PLOIDY_DEFAULT = 1.0;
    double MAX_PLOIDY_DEFAULT = 8;

    static void addOptions(@NotNull Options options) {
        options.addOption(MIN_PURITY, true, "Minimum purity (default " + MIN_PURITY_DEFAULT + ")");
        options.addOption(MAX_PURITY, true, "Maximum purity (default " + MAX_PURITY_DEFAULT + ")");
        options.addOption(PURITY_INCREMENT, true, "Purity increment (default " + PURITY_INCREMENT_DEFAULT + ")");

        options.addOption(MIN_PLOIDY, true, "Minimum ploidy (default " + MIN_PLOIDY_DEFAULT + ")");
        options.addOption(MAX_PLOIDY, true, "Maximum ploidy (default " + MAX_PLOIDY_DEFAULT + ")");
    }

    double minPurity();

    double maxPurity();

    double purityIncrement();

    double minPloidy();

    double maxPloidy();


    @NotNull
    static FittingConfig createConfig(@NotNull final CommandLine cmd) {
        final double minPurity = defaultValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
        final double maxPurity = defaultValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
        final double purityIncrement = defaultValue(cmd, PURITY_INCREMENT, PURITY_INCREMENT_DEFAULT);
        final double minPloidy = defaultValue(cmd, MIN_PLOIDY, MIN_PLOIDY_DEFAULT);
        final double maxPloidy = defaultValue(cmd, MAX_PLOIDY, MAX_PLOIDY_DEFAULT);

        return ImmutableFittingConfig.builder()
                .minPurity(minPurity)
                .maxPurity(maxPurity)
                .purityIncrement(purityIncrement)
                .minPloidy(minPloidy)
                .maxPloidy(maxPloidy)
                .build();

    }

}
