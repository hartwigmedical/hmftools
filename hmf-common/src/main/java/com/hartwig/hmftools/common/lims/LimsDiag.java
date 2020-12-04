package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsDiag {
    COREDB,
    NON_CANCER_DIAGNOSTIC;

    @NotNull
    public static LimsDiag fromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("COREDB")) {
            return COREDB;
        }

        return NON_CANCER_DIAGNOSTIC;
    }
}
