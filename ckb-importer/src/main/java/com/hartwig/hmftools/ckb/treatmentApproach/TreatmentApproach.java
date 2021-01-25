package com.hartwig.hmftools.ckb.treatmentApproach;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApproach {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String profileName();

    @Nullable
    public abstract TreatmentApprochDrugClass drugClass();

    @Nullable
    public abstract TreatmentApproachTherapy therapy();

    @NotNull
    public abstract List<TreatmentApprochReference> reference();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract String updateDate();
}
