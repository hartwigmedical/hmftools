package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsCohortType {
    CORELR02,
    CORERI02,
    CORELR11,
    CORESC11,
    OTHER;

    @NotNull
    public static LimsCohortType fromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("CORELR02")) {
            return LimsCohortType.CORELR02;
        } else if (sampleId.startsWith("CORERI02")) {
            return LimsCohortType.CORERI02;
        } else if (sampleId.startsWith("CORELR11")) {
            return LimsCohortType.CORELR11;
        } else if (sampleId.startsWith("CORESC11")) {
            return LimsCohortType.CORESC11;
        }

        return OTHER;
    }
}
