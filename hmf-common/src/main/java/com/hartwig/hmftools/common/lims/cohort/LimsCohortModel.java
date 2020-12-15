package com.hartwig.hmftools.common.lims.cohort;

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

    @NotNull
    abstract Map<String, LimsCohortConfigData> limsCohortMap();

    private static final Logger LOGGER = LogManager.getLogger(LimsCohortModel.class);

    @Nullable
    public LimsCohortConfigData queryCohortData(@Nullable String cohortString, @NotNull String sampleId) {
        if (cohortString == null) {
            throw new IllegalStateException("Could not resolve LIMS cohort string: '" + cohortString + "'");
        } else {
            LimsCohortConfigData cohortConfigData = limsCohortMap().get(cohortString);
            if (cohortConfigData == null) {
                LOGGER.warn("No cohort map is present for cohortString {}", cohortString);
                return null;
            } else {
                if (sampleId.startsWith(cohortConfigData.cohortId())) {
                    return cohortConfigData;
                } else {
                    LOGGER.error("Cohort '{}' does match with sampleId '{}'", cohortConfigData.cohortId(), sampleId);
                    return null;
                }
            }
        }
    }
}
