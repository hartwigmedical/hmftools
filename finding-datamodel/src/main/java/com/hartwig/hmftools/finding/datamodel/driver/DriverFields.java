package com.hartwig.hmftools.finding.datamodel.driver;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

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
    public boolean isReported()
    {
        return reportedStatus == ReportedStatus.REPORTED;
    }
}
