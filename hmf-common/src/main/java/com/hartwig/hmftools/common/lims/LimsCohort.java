package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum LimsCohort {
    CPCT("CPCT", "CPCT"),
    CPCT_PANCREAS("CPCTPancreas", "CPCT"),
    DRUP("DRUP", "DRUP"),
    DRUP_STAGE3("DRUPstage3", "DRUP"),
    WIDE("WIDE", "WIDE"),
    CORELR02("CORELR02", "CORELR02"),
    CORERI02("CORERI02", "CORERI02"),
    CORELR11("CORELR11", "CORELR11"),
    CORESC11("CORESC11", "CORESC11"),
    COREDB("COREDB", "COREDB"),
    CORE("CORE", "CORE");

    private static final Logger LOGGER = LogManager.getLogger(LimsCohort.class);

    @NotNull
    private final String limsCohortName;
    @NotNull
    private final String expectedSamplePrefix;

    LimsCohort(@NotNull final String limsCohortName, @NotNull final String expectedSamplePrefix) {
        this.limsCohortName = limsCohortName;
        this.expectedSamplePrefix = expectedSamplePrefix;
    }

    @NotNull
    public static LimsCohort fromLimsCohortString(@NotNull String limsCohortString, @NotNull String sampleId) {
        LimsCohort sampleCohort = null;
        for (LimsCohort cohort : LimsCohort.values()) {
            if (cohort.limsCohortName.equals(limsCohortString)) {
                sampleCohort = cohort;
                break;
            }
        }

        if (sampleCohort == null) {
            throw new IllegalStateException("Could not resolve LIMS cohort string: '" + limsCohortString + "'");
        }


        if (!sampleId.startsWith(sampleCohort.expectedSamplePrefix)) {
            LOGGER.error("Cohort '{}' does match with sampleId '{}'", limsCohortString, sampleId);
        }

        return sampleCohort;
    }
}
