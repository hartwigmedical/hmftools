package com.hartwig.hmftools.serve.actionability;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface ActionableEvent {

    @NotNull
    String rawInput();

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
    String urlSource();

    @NotNull
    Set<String> urls();
}
