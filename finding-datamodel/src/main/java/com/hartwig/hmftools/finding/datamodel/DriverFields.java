package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverFields(
        @NotNull String findingKey,
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
