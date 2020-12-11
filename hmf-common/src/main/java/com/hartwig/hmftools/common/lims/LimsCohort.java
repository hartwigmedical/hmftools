package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
    COREDB,
    CORE,
    NON_CANCER;
    private static final Logger LOGGER = LogManager.getLogger(LimsCohort.class);

    @NotNull
    public static LimsCohort fromCohort(@NotNull String cohortString, @NotNull String sampleId) {
        switch (cohortString) {
            case "CPCT":
                if (!sampleId.startsWith("CPCT")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CPCT;
            case "CPCTPancreas":
                if (!sampleId.startsWith("CPCT")) {
                    LOGGER.error("Cohort type and sampleId are don't match");
                }
                return CPCT_PANCREAS;
            case "DRUP":
                if (!sampleId.startsWith("DRUP")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return DRUP;
            case "DRUPstage3":
                if (!sampleId.startsWith("DRUP")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return DRUP_STAGE3;
            case "WIDE":
                if (!sampleId.startsWith("WIDE")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return WIDE;
            case "CORELR02":
                if (!sampleId.startsWith("CORELR02")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CORELR02;
            case "CORERI02":
                if (!sampleId.startsWith("CORERI02")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CORERI02;
            case "CORELR11":
                if (!sampleId.startsWith("CORELR11")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CORELR11;
            case "CORESC11":
                if (!sampleId.startsWith("CORESC11")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CORESC11;
            case "COREDB":
                if (!sampleId.startsWith("COREDB")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return COREDB;
            case "CORE":
                if (!sampleId.startsWith("CORE")) {
                    LOGGER.error("Cohort type {} and sampleId {} are don't match", cohortString, sampleId);
                }
                return CORE;
            default:
                throw new IllegalStateException("Cannot resolve cohort type: " + cohortString);
        }
    }
}
