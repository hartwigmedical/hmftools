package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.cli.Configs.defaultBooleanValue;
import static com.hartwig.hmftools.common.cli.Configs.defaultIntValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BaseQualityRecalibrationConfig {

    String BQR_PLOT = "bqr_plot";
    String BQR_ENABLED = "bqr_enabled";
    String BQR_SAMPLE_SIZE = "bqr_sample_size";
    String BQR_MAX_ALT_COUNT = "bqr_max_alt_count";

    boolean DEFAULT_BQR_PLOT = true;
    boolean DEFAULT_BQR_ENABLED = true;
    int DEFAULT_BQR_MAX_ALT_COUNT = 3;
    int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;

    boolean enabled();

    boolean plot();

    int maxAltCount();

    int sampleSize();

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(BQR_ENABLED, true, "Enabled base quality recalibration [" + DEFAULT_BQR_ENABLED + "]");
        options.addOption(BQR_PLOT, true, "Plot base quality changes  [" + DEFAULT_BQR_PLOT + "]");
        options.addOption(BQR_MAX_ALT_COUNT, true, "Maximum alt count to be an error [" + DEFAULT_BQR_MAX_ALT_COUNT + "]");
        options.addOption(BQR_SAMPLE_SIZE, true, "Sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");
        return options;
    }

    @NotNull
    static BaseQualityRecalibrationConfig createConfig(@NotNull final CommandLine cmd) {
        return ImmutableBaseQualityRecalibrationConfig.builder()
                .enabled(defaultBooleanValue(cmd, BQR_ENABLED, DEFAULT_BQR_ENABLED))
                .plot(defaultBooleanValue(cmd, BQR_PLOT, DEFAULT_BQR_PLOT))
                .maxAltCount(defaultIntValue(cmd, BQR_MAX_ALT_COUNT, DEFAULT_BQR_MAX_ALT_COUNT))
                .sampleSize(defaultIntValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE))
                .build();
    }

}
