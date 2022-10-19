package com.hartwig.hmftools.common.serve.actionability.gene;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.datamodel.gene.GeneLevelEvent;

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
