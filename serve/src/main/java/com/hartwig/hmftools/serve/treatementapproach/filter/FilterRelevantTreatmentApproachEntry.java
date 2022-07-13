package com.hartwig.hmftools.serve.treatementapproach.filter;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FilterRelevantTreatmentApproachEntry {

    @NotNull
    public abstract String treatment();

    @NotNull
    public abstract String eventMatch();

    @Nullable
    public abstract EvidenceLevel level();

    @NotNull
    public abstract EvidenceDirection direction();
}
