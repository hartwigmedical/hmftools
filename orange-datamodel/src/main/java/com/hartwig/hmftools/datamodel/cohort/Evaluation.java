package com.hartwig.hmftools.datamodel.cohort;

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
