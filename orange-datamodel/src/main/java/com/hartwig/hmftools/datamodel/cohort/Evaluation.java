package com.hartwig.hmftools.datamodel.cohort;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Evaluation
{
    @Nullable
    String cancerType();

    double panCancerPercentile();

    @Nullable
    Double cancerTypePercentile();
}
