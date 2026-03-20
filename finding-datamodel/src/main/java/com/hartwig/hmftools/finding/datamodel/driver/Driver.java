package com.hartwig.hmftools.finding.datamodel.driver;

import java.util.Set;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.ReportedStatus;

import jakarta.validation.constraints.NotNull;

public interface Driver extends Finding
{
    @NotNull DriverSource driverSource();

    @NotNull
    ReportedStatus reportedStatus();

    @NotNull
    DriverInterpretation driverInterpretation();

    double driverLikelihood();

    boolean isReported();

    Set<String> genes();
}
