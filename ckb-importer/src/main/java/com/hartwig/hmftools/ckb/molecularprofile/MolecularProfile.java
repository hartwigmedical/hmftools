package com.hartwig.hmftools.ckb.molecularprofile;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfile {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<MolecularProfileGeneVariant> geneVariant();

    @NotNull
    public abstract List<MolecularProfileProfileTreatmentApproache> profileProfileTreatmentApproache();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract String updateDate();

    @NotNull
    public abstract MolecularProfileComplexMolecularProfileEvidence complexMolecularProfileEvidence();

}
