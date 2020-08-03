package com.hartwig.hmftools.serve.vicc.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FusionAnnotation {

    @NotNull
    public abstract String fusion();

    @NotNull
    public abstract String fusionName();
}
