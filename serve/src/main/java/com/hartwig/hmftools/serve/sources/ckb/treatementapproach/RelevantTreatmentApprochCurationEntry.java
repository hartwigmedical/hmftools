package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RelevantTreatmentApprochCurationEntry {

    @NotNull
    public abstract RelevantTreatmentApproachCurationType curationType();

    @NotNull
    public abstract RelevantTreatmentApprochCurationEntryKey curationKey();

    @NotNull
    public abstract String curatedtreatmentApproach();
}