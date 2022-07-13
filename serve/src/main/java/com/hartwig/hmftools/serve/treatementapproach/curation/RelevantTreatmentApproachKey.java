package com.hartwig.hmftools.serve.treatementapproach.curation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RelevantTreatmentApproachKey {

    @NotNull
    public abstract String treatment();

    @NotNull
    public abstract String treatmentApproach();

    @NotNull
    public abstract String matchEvent();

}
