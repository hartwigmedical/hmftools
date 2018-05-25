package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.jetbrains.annotations.NotNull;

public final class CopiesOnReportDependentOnPurity {

    private CopiesOnReportDependentOnPurity() {
    }

    @NotNull
    public static String copiesValue(@NotNull final FittedPurityStatus fitStatus, @NotNull final String copiesValue) {
        return fitStatus == FittedPurityStatus.NO_TUMOR ? "N/A" : copiesValue;
    }
}
