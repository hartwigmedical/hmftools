package com.hartwig.hmftools.finding.datamodel.driver;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.ReportedStatus;

import jakarta.validation.constraints.NotNull;

public interface Driver extends Finding
{
    @NotNull DriverSource driverSource();

    @NotNull
    ReportedStatus reportedStatus();

    @NotNull DriverInterpretation driverInterpretation();

    double driverLikelihood();

    boolean isReported();

    void accept(@NotNull DriverVisitor visitor);
}
