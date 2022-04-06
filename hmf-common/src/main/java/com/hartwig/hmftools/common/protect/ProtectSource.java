package com.hartwig.hmftools.common.protect;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ProtectSource {

    @NotNull
    public abstract Knowledgebase sources();

    @NotNull
    public abstract String sourceEvent();

    @NotNull
    public abstract Set<String> sourceUrls();

    @NotNull
    public abstract ProtectEvidenceType evidenceType();

    @Nullable
    public abstract Integer rangeRank();
}
