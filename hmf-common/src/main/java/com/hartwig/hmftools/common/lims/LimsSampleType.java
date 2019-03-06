package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsSampleType {
    CORE,
    WIDE,
    CPCT,
    DRUP,
    COLO;

    @NotNull
    public static LimsSampleType fromSampleId(@NotNull String sampleId) {
        // TODO
        return CORE;
    }
}
