package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsSampleType {
    CORE,
    WIDE,
    CPCT,
    DRUP,
    COLO,
    PNT;

    @NotNull
    public static LimsSampleType fromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("CPCT")) {
            return CPCT;
        } else if (sampleId.startsWith("WIDE")) {
            return WIDE;
        } else if (sampleId.startsWith("CORE")) {
            return CORE;
        } else if (sampleId.startsWith("DRUP")) {
            return DRUP;
        } else if (sampleId.startsWith("COLO")) {
            return COLO;
        } else if (sampleId.startsWith("PNT")) { // is only used for test rapport COLO for extern uses
            return PNT;
        }

        throw new IllegalStateException("Cannot resolve type for sampleId: " + sampleId);
    }
}
