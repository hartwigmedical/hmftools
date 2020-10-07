package com.hartwig.hmftools.serve.actionability;

import com.hartwig.hmftools.serve.Source;

import org.jetbrains.annotations.NotNull;

public abstract class ActionableEvent {

    @NotNull
    public abstract Source source();

    @NotNull
    public abstract String treatment();

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
