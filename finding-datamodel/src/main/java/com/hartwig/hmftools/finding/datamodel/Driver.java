package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

public interface Driver extends Finding
{
    @NotNull DriverSource driverSource();

    @NotNull ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

    double driverLikelihood();
}
