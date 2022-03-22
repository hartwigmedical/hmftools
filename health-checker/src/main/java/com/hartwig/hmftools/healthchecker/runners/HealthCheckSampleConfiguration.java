package com.hartwig.hmftools.healthchecker.runners;

import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
public interface HealthCheckSampleConfiguration {

    String sampleName();

    String wgsMetricsFile();

    String flagstatFile();

    static HealthCheckSampleConfiguration of(@Nullable String sampleName, @Nullable String wgsMetricsFile, @Nullable String flagstatFile) {
        if (sampleName != null) {
            return ImmutableHealthCheckSampleConfiguration.builder()
                    .sampleName(sampleName)
                    .wgsMetricsFile(notnullOrThrow(wgsMetricsFile))
                    .flagstatFile(notnullOrThrow(flagstatFile))
                    .build();
        }
        return null;
    }

    @NotNull
    static String notnullOrThrow(@Nullable final String wgsMetricsFile) {
        return Optional.ofNullable(wgsMetricsFile).orElseThrow();
    }
}
