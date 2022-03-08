package com.hartwig.hmftools.serve.actionability;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklisting;

import org.jetbrains.annotations.NotNull;

public interface ActionableEvent2 {

    @NotNull
    Knowledgebase source();

    @NotNull
    String sourceEvent();

    @NotNull
    Set<String> sourceUrls();

    @NotNull
    String treatment();

    @NotNull
    String cancerType();

    @NotNull
    String doid();

    @NotNull
    Set<TumorLocationBlacklisting> blacklistings();

    @NotNull
    EvidenceLevel level();

    @NotNull
    EvidenceDirection direction();

    @NotNull
    Set<String> evidenceUrls();
}
