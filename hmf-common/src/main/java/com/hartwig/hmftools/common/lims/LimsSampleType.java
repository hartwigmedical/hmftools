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
        }

        throw new IllegalStateException("Cannot resolve type for sampleId: " + sampleId);
    }

    @NotNull
    public static LimsSampleType fromRunName(@NotNull String runName) {
        if (runName.contains("CPCT")) {
            return CPCT;
        } else if (runName.contains("WIDE")) {
            return WIDE;
        } else if (runName.contains("CORE")) {
            return CORE;
        } else if (runName.contains("DRUP")) {
            return DRUP;
        } else if (runName.contains("COLO")) {
            return COLO;
        }

        throw new IllegalStateException("Cannot resolve type for run name: " + runName);
    }
}
