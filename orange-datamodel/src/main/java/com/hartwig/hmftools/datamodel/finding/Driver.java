package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

public interface Driver extends Finding
{
    boolean isReportable();
    boolean isCandidate();
    @NotNull
    DriverInterpretation driverInterpretation();
}
