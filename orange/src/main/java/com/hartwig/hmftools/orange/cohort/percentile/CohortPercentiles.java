package com.hartwig.hmftools.orange.cohort.percentile;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CohortPercentiles {

    @NotNull
    public abstract String cancerType();

    @NotNull 
    public abstract List<Double> percentileValues();

}
