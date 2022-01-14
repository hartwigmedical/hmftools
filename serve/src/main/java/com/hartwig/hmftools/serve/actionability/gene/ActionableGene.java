package com.hartwig.hmftools.serve.actionability.gene;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableGene implements ActionableEvent {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract GeneLevelEvent event();

}
