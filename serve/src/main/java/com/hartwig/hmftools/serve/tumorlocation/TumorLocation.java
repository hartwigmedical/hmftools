package com.hartwig.hmftools.serve.tumorlocation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorLocation {

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String doid();
}
