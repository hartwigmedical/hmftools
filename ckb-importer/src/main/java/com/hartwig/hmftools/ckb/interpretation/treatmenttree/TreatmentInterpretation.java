package com.hartwig.hmftools.ckb.interpretation.treatmenttree;

import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentInterpretation {

    @NotNull
    public abstract TreatmentApprochInterpretation treatmentApproach();


}
