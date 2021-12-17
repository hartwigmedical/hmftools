package com.hartwig.hmftools.serve.sources.actin.filter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinFilterEntry {

    @NotNull
    public abstract ActinFilterType type();

    @NotNull
    public abstract String value();
}
