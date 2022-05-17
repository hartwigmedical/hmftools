package com.hartwig.hmftools.summon.actionability;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityEntry {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract Type type();

    @NotNull
    public abstract String conclusion();
}