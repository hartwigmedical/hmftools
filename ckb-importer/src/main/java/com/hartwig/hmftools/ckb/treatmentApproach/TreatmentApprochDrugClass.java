package com.hartwig.hmftools.ckb.treatmentApproach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApprochDrugClass {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String drugClass();
}
