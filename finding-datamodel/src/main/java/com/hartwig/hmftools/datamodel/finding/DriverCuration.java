package com.hartwig.hmftools.datamodel.finding;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverCuration(
        @NotNull String findingKey,
        long utcTimestamp,
        @NotNull String userId,
        boolean reportableOverride,
        @NotNull String comment)
{
}
