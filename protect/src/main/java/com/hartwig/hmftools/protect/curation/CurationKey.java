package com.hartwig.hmftools.protect.curation;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class CurationKey {

    @NotNull
    public abstract String eventKeyword();

    @NotNull
    public abstract String treatment();

    @NotNull
    public abstract EvidenceDirection direction();

}
