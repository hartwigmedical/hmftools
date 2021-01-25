package com.hartwig.hmftools.ckb.variant;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantMolecularProfile {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String profileName();

    @Nullable
    public abstract List<VariantProfileTreatmentApproach> profileTreatmentApproach();
}
