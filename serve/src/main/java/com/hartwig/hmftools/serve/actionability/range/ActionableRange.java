package com.hartwig.hmftools.serve.actionability.range;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableRange extends ActionableEvent {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String end();

    @NotNull
    public abstract String mutationType();
}
