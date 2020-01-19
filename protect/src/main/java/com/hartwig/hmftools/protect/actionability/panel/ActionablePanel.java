package com.hartwig.hmftools.protect.actionability.panel;

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

    public abstract boolean variant();

    public abstract boolean drup();

    @NotNull
    public abstract String responsive();

    @NotNull
    public abstract String responsiveSource();

    @NotNull
    public abstract String resistant();

    @NotNull
    public abstract String resistantSource();

    @Override
    public int compareTo(@NotNull final ActionablePanel o) {
        return gene().compareTo(o.gene());
    }
}
