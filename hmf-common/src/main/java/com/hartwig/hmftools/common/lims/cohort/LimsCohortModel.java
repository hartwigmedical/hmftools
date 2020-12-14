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
    abstract Map<String, LimsCohortConfigData> limsCohortMap();

    @Nullable
    public LimsCohortConfigData queryCohortData(@Nullable String cohortString) {
        if (cohortString == null) {
            return null;
        }
        LimsCohortConfigData configData = limsCohortMap().get(cohortString);
        return configData;


    }


}
