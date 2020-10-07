package com.hartwig.hmftools.serve.actionability;

import org.jetbrains.annotations.NotNull;

public abstract class ActionableEvent {

    @NotNull
    public abstract String source();

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
}
