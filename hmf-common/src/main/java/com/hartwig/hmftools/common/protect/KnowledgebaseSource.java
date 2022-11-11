package com.hartwig.hmftools.common.protect;

import java.util.Set;

import com.hartwig.serve.datamodel.Knowledgebase;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnowledgebaseSource {

    @NotNull
    public abstract Knowledgebase name();

    @NotNull
    public abstract String sourceEvent();

    @NotNull
    public abstract Set<String> sourceUrls();

    @NotNull
    public abstract EvidenceType evidenceType();

    @Nullable
    public abstract Integer rangeRank();

    @NotNull
    public abstract Set<String> evidenceUrls();
}
