package com.hartwig.hmftools.serve.actionability.fusion;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableFusion implements ActionableEvent {

    @NotNull
    public abstract String geneUp();

    @Nullable
    public abstract Integer exonUp();

    @NotNull
    public abstract String geneDown();

    @Nullable
    public abstract Integer exonDown();
}