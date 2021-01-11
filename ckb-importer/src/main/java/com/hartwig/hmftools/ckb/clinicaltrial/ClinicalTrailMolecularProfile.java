package com.hartwig.hmftools.ckb.clinicaltrial;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrailMolecularProfile {

    @Nullable
    public abstract String id();

    @NotNull
    public abstract String profileName();
}
