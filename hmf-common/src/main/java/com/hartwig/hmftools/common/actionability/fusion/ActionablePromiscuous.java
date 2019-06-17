package com.hartwig.hmftools.common.actionability.fusion;

import com.hartwig.hmftools.common.actionability.Actionable;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionablePromiscuous implements Actionable {
    @NotNull
    public abstract String gene();
}
