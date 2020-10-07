package com.hartwig.hmftools.serve.vicc;

import com.hartwig.hmftools.serve.actionability.EvidenceDirection;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableEvidence {

    @NotNull
    public abstract String drugs();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String doid();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract EvidenceDirection direction();

    @NotNull
    public abstract String url();
}
