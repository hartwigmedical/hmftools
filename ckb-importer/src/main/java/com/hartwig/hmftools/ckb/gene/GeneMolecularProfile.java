package com.hartwig.hmftools.ckb.gene;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneMolecularProfile {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<GeneProfileTreatmentApproache> profileTreatmentApproache();
}
