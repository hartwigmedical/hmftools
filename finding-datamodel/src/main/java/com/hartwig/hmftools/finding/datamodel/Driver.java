package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;

import jakarta.validation.constraints.NotNull;

public interface Driver extends Event
{
    @NotNull DriverSource driverSource();

    @NotNull ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

    double driverLikelihood();

    boolean isReported();

    Set<String> genes();
}
