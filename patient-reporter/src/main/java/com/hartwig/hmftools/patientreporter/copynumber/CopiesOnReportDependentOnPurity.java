package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.jetbrains.annotations.NotNull;

public abstract class CopiesOnReportDependentOnPurity {

    public PurpleAnalysis purple;

    @NotNull
    public String copiesValue(final int value) {
        return purple.status() == FittedPurityStatus.NO_TUMOR ? "N/A" : Integer.toString(value);
    }
}
