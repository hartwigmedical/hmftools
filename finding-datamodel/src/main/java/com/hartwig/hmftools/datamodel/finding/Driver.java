package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;

import jakarta.validation.constraints.NotNull;

public interface Driver extends Finding
{
    @NotNull DriverSource driverSource();

    @NotNull ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

    double driverLikelihood();
}
