package com.hartwig.hmftools.serve.actionability;

import java.util.Set;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.Sources;

import org.jetbrains.annotations.NotNull;

public interface ActionableEvent {

    @NotNull
    Sources source();

    @NotNull
    String treatment();

    @NotNull
    String cancerType();

    @NotNull
    String doid();

    @NotNull
    String tumorLocationBlacklisting();

    @NotNull
    EvidenceLevel level();

    @NotNull
    EvidenceDirection direction();

    @NotNull
    Set<String> sourceUrls();

    @NotNull
    Set<String> evidenceUrls();
}