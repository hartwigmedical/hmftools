package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RelevantTreatmentApprochCurationEntryKey {

    @NotNull
    public abstract String treatment();

    @Nullable
    public abstract String treatmentApproach();

    @NotNull
    public abstract String event();

    @NotNull
    public abstract EvidenceDirection direction();
}