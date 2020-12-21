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
    public LimsCohortConfig queryCohortData(@Nullable String cohortString, @NotNull String sampleId, @NotNull String sampleIdMetadata) {
        String cohortStringNew = !sampleIdMetadata.startsWith("COLO") ? cohortString : "COLO";

        if (cohortStringNew == null) {
            throw new IllegalStateException("No cohort string present in LIMS for sample '{}'" + sampleId);
        } else {
            LimsCohortConfig cohortConfigData = limsCohortMap().get(cohortStringNew);
            if (cohortConfigData == null) {
                throw new IllegalStateException("No cohort config present for cohort '{}'" + cohortStringNew);
            } else {
                if (cohortConfigData.cohortId().contains("CPCT") || cohortConfigData.cohortId().contains("DRUP")) {
                    if (sampleId.startsWith(cohortConfigData.cohortId().substring(0, 4))
                            || sampleIdMetadata.startsWith(cohortConfigData.cohortId().substring(0, 4))) {
                        return cohortConfigData;
                    } else {
                        throw new IllegalStateException("No matching for DRUP/CPCT");
                    }
                } else if (!cohortConfigData.cohortId().contains("CPCT") || !cohortConfigData.cohortId().contains("DRUP")) {
                    return cohortConfigData;
                } else {
                    throw new IllegalStateException(
                            "Cohort '{}' does not seem to match with sample '{}'" + cohortConfigData.cohortId() + sampleId);
                }
            }
        }
    }
}
