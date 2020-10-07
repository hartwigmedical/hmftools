package com.hartwig.hmftools.serve.actionability.hotspot;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableHotspot extends ActionableEvent {

    @NotNull
    public abstract String chromosome();

    public abstract long position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();
}


