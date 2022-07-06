package com.hartwig.hmftools.serve.curation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RelevantTreatmentApproch {

    @NotNull
    public abstract RelevantTreatmentApproachKey treatmentApproachKey();

    @NotNull
    public abstract String curatedtreatmentApproach();
}