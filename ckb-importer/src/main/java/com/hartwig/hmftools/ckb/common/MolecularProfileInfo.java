package com.hartwig.hmftools.ckb.common;

import java.util.List;

import com.hartwig.hmftools.ckb.gene.GeneProfileTreatmentApproache;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileInfo {

    @Nullable
    public abstract String id();

    @NotNull
    public abstract String profileName();

    public abstract List<GeneProfileTreatmentApproache> profileTreatmentApproache();

}
