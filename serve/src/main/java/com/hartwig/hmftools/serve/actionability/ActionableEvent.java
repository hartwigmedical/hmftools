package com.hartwig.hmftools.serve.actionability;

import com.hartwig.hmftools.common.serve.EvidenceDirection;
import com.hartwig.hmftools.common.serve.EvidenceLevel;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public interface ActionableEvent {

    @NotNull
    Knowledgebase source();

    @NotNull
    String treatment();

    @NotNull
    String cancerType();

    @NotNull
    String doid();

    @NotNull
    EvidenceLevel level();

    @NotNull
    EvidenceDirection direction();

    @NotNull
    String url();
}
