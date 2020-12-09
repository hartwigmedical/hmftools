package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsCohort {
    CPCT,
    CPCT_PANCREAS,
    DRUP,
    DRUP_STAGE3,
    WIDE,
    CORELR02,
    CORERI02,
    CORELR11,
    CORESC11,
    COREDB01,
    CORE,
    NON_CANCER;

    @NotNull
    public static LimsCohort fromCohort(@NotNull String cohortString) {
        switch (cohortString) {
            case "CPCT":
                return CPCT;
            case "CPCTpancreas":
                return CPCT_PANCREAS;
            case "DRUP":
                return DRUP;
            case "DRUPstage3":
                return DRUP_STAGE3;
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
                throw new IllegalStateException("Cannot resolve cohort type: " + cohortString);
        }
    }
}
