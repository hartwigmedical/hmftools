package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsCohort {
    CPCT,
    DRUP,
    WIDE,
    CORELR02,
    CORERI02,
    CORELR11,
    CORESC11,
    COREDB01,
    CORE,
    NON_CANCER;

    @NotNull
    public static LimsCohort fromCohort(@NotNull String cohortString, @NotNull String sampleId) {
        switch (cohortString) {
            case "CPCT":
            case "CPCTpancreas":
                return CPCT;
            case "DRUP":
            case "DRUPstage3":
                return DRUP;
            case "WIDE":
                return WIDE;
            case "CORE":
                return CORE;
            case "CORELR02":
                return CORELR02;
            case "CORERI02":
                return CORERI02;
            case "CORELR11":
                return CORELR11;
            case "CORESC11":
                return CORESC11;
            case "COREDB01":
                return COREDB01;

            default:
                throw new IllegalStateException(
                        "Cannot resolve cohort type: " + cohortString+  " of sample: " + sampleId );
        }

    }

}
