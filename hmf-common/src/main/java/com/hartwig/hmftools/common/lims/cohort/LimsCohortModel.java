package com.hartwig.hmftools.common.lims.cohort;

import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsCohortModel {

    @NotNull
    protected abstract Map<String, LimsCohortConfig> limsCohortMap();

    @Nullable
    public LimsCohortConfig queryCohortData(@Nullable String cohortString, @NotNull String sampleId) {
        if (cohortString == null) {
            throw new IllegalStateException("No cohort string present in LIMS for sample '{}'" + sampleId);
        } else {
            LimsCohortConfig cohortConfigData = limsCohortMap().get(cohortString);
            if (cohortConfigData == null) {
                throw new IllegalStateException("No cohort config present for cohort '{}'" + cohortString);
            } else {
                if (sampleId.startsWith(cohortConfigData.cohortId())) {
                    return cohortConfigData;
                } else {
                    throw new IllegalStateException(
                            "Cohort '{}' does not seem to match with sample '{}'" + cohortConfigData.cohortId() + sampleId);
                }
            }
        }
    }
}
