package com.hartwig.hmftools.serve.actionability.immuno;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableHLA implements ActionableEvent {

    @NotNull
    public abstract String hlaTypering();
}
