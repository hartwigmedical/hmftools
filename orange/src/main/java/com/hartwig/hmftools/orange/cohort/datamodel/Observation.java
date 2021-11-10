package com.hartwig.hmftools.orange.cohort.datamodel;

import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Observation {

    @NotNull
    public abstract PercentileType type();

    @NotNull
    public abstract Sample sample();

    public abstract double value();
}
