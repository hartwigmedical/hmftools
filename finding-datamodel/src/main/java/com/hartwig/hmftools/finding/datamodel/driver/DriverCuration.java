package com.hartwig.hmftools.finding.datamodel.driver;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

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
