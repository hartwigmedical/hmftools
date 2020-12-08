package com.hartwig.hmftools.serve.extraction.gene;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneLevelAnnotation {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract GeneLevelEvent event();
}
