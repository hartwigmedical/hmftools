package com.hartwig.hmftools.common.actionability.panel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionablePanel implements Comparable<ActionablePanel> {

    @NotNull
    public abstract String gene();

    public abstract boolean amplification();

    public abstract boolean deletion();

    public abstract boolean fusion();

    public abstract boolean variant();

    @Override
    public int compareTo(@NotNull final ActionablePanel o) {
        return gene().compareTo(o.gene());
    }
}
