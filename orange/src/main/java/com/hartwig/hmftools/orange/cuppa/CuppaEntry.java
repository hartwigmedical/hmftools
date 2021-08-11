package com.hartwig.hmftools.orange.cuppa;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuppaEntry {

    @NotNull
    public abstract String category();

    @NotNull
    public abstract String resultType();

    @NotNull
    public abstract String dataType();

    @NotNull
    public abstract String value();

    @NotNull
    public abstract String refCancerType();

    public abstract double refValue();
}
