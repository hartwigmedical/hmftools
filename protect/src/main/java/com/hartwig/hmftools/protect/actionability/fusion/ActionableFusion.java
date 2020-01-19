package com.hartwig.hmftools.protect.actionability.fusion;

import com.hartwig.hmftools.protect.actionability.Actionable;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableFusion implements Actionable {

    @NotNull
    public abstract String fiveGene();

    @NotNull
    public abstract String threeGene();

}
