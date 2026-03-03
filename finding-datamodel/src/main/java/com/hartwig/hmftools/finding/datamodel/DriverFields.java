package com.hartwig.hmftools.finding.datamodel;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverFields(
        @NotNull String findingKey,
        @NotNull String event,
        @NotNull DriverSource driverSource,
        @NotNull ReportedStatus reportedStatus,
        @NotNull DriverInterpretation driverInterpretation,
        double driverLikelihood
)
{
    boolean isReported()
    {
        return reportedStatus == ReportedStatus.REPORTED;
    }
}
