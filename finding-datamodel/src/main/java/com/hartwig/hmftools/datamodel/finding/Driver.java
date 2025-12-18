package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.jetbrains.annotations.NotNull;

public interface Driver extends Finding
{
    @NotNull DriverSource driverSource();

    @NotNull ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

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
