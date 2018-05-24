package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;

import org.jetbrains.annotations.NotNull;

public final class CopiesOnReportDependentOnPurity {

    static PurpleAnalysis purple;

    @NotNull
    public static String copiesValue(final String value) {
        return purple.status() == FittedPurityStatus.NO_TUMOR ? "N/A" : value;
    }
}
