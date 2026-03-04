package com.hartwig.hmftools.finding.datamodel;

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
