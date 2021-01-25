package com.hartwig.hmftools.ckb.molecularprofile;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfileExtendedEvidenceList {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract MolecularProfileMolecularProfile molecularProfile();

    @NotNull
    public abstract MolecularProfileTherapy therapy();

    @NotNull
    public abstract MolecularProfileIndication indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<MolecularProfileReferences> reference();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();
}
