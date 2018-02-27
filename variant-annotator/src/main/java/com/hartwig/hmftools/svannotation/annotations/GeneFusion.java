package com.hartwig.hmftools.svannotation.annotations;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneFusion {

    public abstract boolean reportable();

    @NotNull
    public abstract Transcript upstreamLinkedAnnotation();

    @NotNull
    public abstract Transcript downstreamLinkedAnnotation();

    @NotNull
    public abstract String cosmicURL();
}