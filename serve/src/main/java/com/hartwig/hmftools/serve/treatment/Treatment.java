package com.hartwig.hmftools.serve.treatment;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Treatment {

    @NotNull
    public abstract String treament();

    @NotNull
    public abstract Set<String> drugClasses();
}
