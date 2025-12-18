package com.hartwig.hmftools.datamodel.driver;

import com.hartwig.hmftools.datamodel.finding.Finding;

import org.jetbrains.annotations.NotNull;

public interface Driver extends Finding
{
    @NotNull
    DriverSource driverSource();

    @NotNull
    ReportedStatus reportedStatus();

    @NotNull
    DriverInterpretation driverInterpretation();

    default boolean isReportable()
    {
        return reportedStatus() == ReportedStatus.REPORTED &&
                (driverInterpretation() == DriverInterpretation.HIGH || driverInterpretation() == DriverInterpretation.MEDIUM);
    }

    default boolean isCandidate()
    {
        return reportedStatus() == ReportedStatus.REPORTED && driverInterpretation() == DriverInterpretation.LOW;
    }
}
