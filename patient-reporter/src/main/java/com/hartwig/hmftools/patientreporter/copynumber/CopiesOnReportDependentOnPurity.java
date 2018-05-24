package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.jetbrains.annotations.NotNull;

public final class CopiesOnReportDependentOnPurity {

    @NotNull
    public static String copiesValue(@NotNull final String fitStatus, @NotNull final String valueOfTotalCopies) {
        return fitStatus.equals(FittedPurityStatus.NO_TUMOR) ? "N/A" : valueOfTotalCopies;
    }
}
