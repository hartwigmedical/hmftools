package com.hartwig.hmftools.serve.actionability;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.tumorlocation.TumorLocation;

import org.jetbrains.annotations.NotNull;

public interface ActionableEvent {

    @NotNull
    Knowledgebase source();

    @NotNull
    String sourceEvent();

    @NotNull
    Set<String> sourceUrls();

    @NotNull
    String treatment();

    @NotNull
    TumorLocation whiteListCancerType();

    @NotNull
    Set<TumorLocation> blackListCancerTypes();

    @NotNull
    EvidenceLevel level();

    @NotNull
    EvidenceDirection direction();

    @NotNull
    Set<String> evidenceUrls();
}