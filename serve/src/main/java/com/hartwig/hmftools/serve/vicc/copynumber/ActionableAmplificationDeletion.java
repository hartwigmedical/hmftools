package com.hartwig.hmftools.serve.vicc.copynumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableAmplificationDeletion {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String eventType();

    @NotNull
    public abstract String source();

    //Below items only present when it is a actionable amplification or deletion
    @Nullable
    public abstract String drug();

    @Nullable
    public abstract String drugType();

    @Nullable
    public abstract String cancerType();

    @Nullable
    public abstract String level();

    @Nullable
    public abstract String direction();

    @Nullable
    public abstract String sourceLink();
}