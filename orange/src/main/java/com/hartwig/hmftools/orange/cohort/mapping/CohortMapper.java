package com.hartwig.hmftools.orange.cohort.mapping;

import com.hartwig.hmftools.orange.cohort.datamodel.Sample;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface CohortMapper {

    @Nullable
    String cancerTypeForSample(@NotNull Sample sample);
}
