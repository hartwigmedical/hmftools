package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsCoreCohort {
    CORELR02,
    CORERI02,
    CORELR11,
    CORESC11,
    OTHER;

    @NotNull
    public static LimsCoreCohort fromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("CORELR02")) {
            return LimsCoreCohort.CORELR02;
        } else if (sampleId.startsWith("CORERI02")) {
            return LimsCoreCohort.CORERI02;
        } else if (sampleId.startsWith("CORELR11")) {
            return LimsCoreCohort.CORELR11;
        } else if (sampleId.startsWith("CORESC11")) {
            return LimsCoreCohort.CORESC11;
        }

        return OTHER;
    }
}
