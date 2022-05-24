package com.hartwig.hmftools.rose.actionability;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityKey {

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract TypeAlteration type();
}
