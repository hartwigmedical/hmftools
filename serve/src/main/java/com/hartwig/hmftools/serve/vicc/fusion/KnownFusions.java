package com.hartwig.hmftools.serve.vicc.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownFusions {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String eventType();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String sourceLink();
}

