package com.hartwig.hmftools.datamodel.driver;

import com.hartwig.hmftools.datamodel.finding.Finding;

import org.immutables.gson.Gson;
import org.jetbrains.annotations.NotNull;

public interface Driver extends Finding
{
    @NotNull ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

    @Gson.Ignore
    default boolean isReported()
    {
        return reportedStatus() == ReportedStatus.REPORTED &&
                (driverInterpretation() == DriverInterpretation.HIGH || driverInterpretation() == DriverInterpretation.MEDIUM);
    }

    @Gson.Ignore
    default boolean isCandidate()
    {
        return reportedStatus() == ReportedStatus.CANDIDATE &&
                driverInterpretation() == DriverInterpretation.LOW;
    }
}
