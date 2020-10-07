package com.hartwig.hmftools.serve.actionability;

import com.hartwig.hmftools.serve.Source;

import org.jetbrains.annotations.NotNull;

public interface ActionableEvent {

    @NotNull
    Source source();

    @NotNull
    String treatment();

    @NotNull
    String cancerType();

    @NotNull
    String doid();

    @NotNull
    String level();

    @NotNull
    EvidenceDirection direction();

    @NotNull
    String url();
}
