package com.hartwig.hmftools.orange.cohort.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Evaluation {

    @Nullable
    public abstract String cancerType();

    public abstract double panCancerPercentile();

    @Nullable
    public abstract Double cancerTypePercentile();
}
