package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsSampleType {
    CORE,
    WIDE,
    CPCT,
    DRUP,
    OTHER;

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
        } else if (sampleId.startsWith("COLO") || sampleId.startsWith("PNT")) {
            // PNT is only used for test rapport COLO for extern uses
            return OTHER;
        }

        throw new IllegalStateException("Cannot resolve type for sampleId: " + sampleId);
    }
}
