package com.hartwig.hmftools.patientdb.clinical.lims.cohort;

import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsCohortModel {

    private static final Logger LOGGER = LogManager.getLogger(LimsCohortModel.class);

    @NotNull
    protected abstract Map<String, LimsCohortConfig> limsCohortMap();

    @Nullable
    public LimsCohortConfig queryCohortData(@Nullable String cohortString, @NotNull String sampleId) {
        if (cohortString == null || cohortString.isEmpty()) {
            LOGGER.warn("No cohort string present in LIMS for sample '{}'.", sampleId);
            return null;
        }
        LimsCohortConfig cohortConfigData = limsCohortMap().get(cohortString);
        if (cohortConfigData == null) {
            LOGGER.warn("Could not resolve cohort config for sample '{}' based on LIMS cohort '{}'", sampleId, cohortString);
            return null;
        }
        return cohortConfigData;
    }
}