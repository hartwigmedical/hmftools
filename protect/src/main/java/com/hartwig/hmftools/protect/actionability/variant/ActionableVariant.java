package com.hartwig.hmftools.protect.actionability.variant;

import com.hartwig.hmftools.protect.actionability.Actionable;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableVariant implements Actionable {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    public abstract long position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

}
