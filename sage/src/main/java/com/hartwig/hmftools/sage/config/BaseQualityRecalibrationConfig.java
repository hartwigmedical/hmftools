package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

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
    String BQR_MIN_MAP_QUAL = "bqr_min_map_qual";

    boolean DEFAULT_BQR_PLOT = true;
    boolean DEFAULT_BQR_ENABLED = true;
    int DEFAULT_BQR_MAX_ALT_COUNT = 3;
    int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;
    int DEFAULT_BQR_MIN_MAP_QUAL = 10;

    boolean enabled();

    boolean plot();

    int maxAltCount();

    int sampleSize();

    int minMapQuality();

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(BQR_ENABLED, true, "BQR (Base Quality Recalibration) enabled [" + DEFAULT_BQR_ENABLED + "]");
        options.addOption(BQR_PLOT, true, "BQR plots [" + DEFAULT_BQR_PLOT + "]");
        options.addOption(BQR_MAX_ALT_COUNT, true, "BQR maximum alt count to be an error [" + DEFAULT_BQR_MAX_ALT_COUNT + "]");
        options.addOption(BQR_SAMPLE_SIZE, true, "BQR sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");
        options.addOption(BQR_MIN_MAP_QUAL, true, "BQR min base quality remap qual [" + DEFAULT_BQR_MIN_MAP_QUAL + "]");
        return options;
    }

    @NotNull
    static BaseQualityRecalibrationConfig createConfig(@NotNull final CommandLine cmd) {
        return ImmutableBaseQualityRecalibrationConfig.builder()
                .enabled(getConfigValue(cmd, BQR_ENABLED, DEFAULT_BQR_ENABLED))
                .plot(getConfigValue(cmd, BQR_PLOT, DEFAULT_BQR_PLOT))
                .maxAltCount(getConfigValue(cmd, BQR_MAX_ALT_COUNT, DEFAULT_BQR_MAX_ALT_COUNT))
                .sampleSize(getConfigValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE))
                .minMapQuality(getConfigValue(cmd, BQR_MIN_MAP_QUAL, DEFAULT_BQR_MIN_MAP_QUAL))
                .build();
    }
}
