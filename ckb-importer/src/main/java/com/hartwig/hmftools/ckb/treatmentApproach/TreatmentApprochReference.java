package com.hartwig.hmftools.ckb.treatmentApproach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApprochReference {

    @NotNull
    public abstract String id();

    @Nullable
    public abstract String pubMedId();

    @NotNull
    public abstract String title();

    @Nullable
    public abstract String url();
}
