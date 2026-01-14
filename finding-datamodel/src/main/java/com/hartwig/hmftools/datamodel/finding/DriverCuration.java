package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

public record DriverCuration(
        @NotNull String findingKey,
        long utcTimestamp,
        @NotNull String userId,
        boolean reportableOverride,
        @NotNull String comment) {}
