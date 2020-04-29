package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsStudy {
    CORE,
    WIDE,
    CPCT,
    DRUP,
    NON_STUDY;

    @NotNull
    public static LimsStudy fromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("CPCT")) {
            return CPCT;
        } else if (sampleId.startsWith("WIDE")) {
            return WIDE;
        } else if (sampleId.startsWith("CORE")) {
            return CORE;
        } else if (sampleId.startsWith("DRUP")) {
            return DRUP;
        }

        return NON_STUDY;
    }
}
