package com.hartwig.hmftools.ckb.datamodelinterpretation.molecularprofile;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularProfile {

    public abstract int id();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract Date createDate();

    @NotNull
    public abstract Date updateDate();

    @NotNull
    public abstract List<MolecularProfileClinicalTrial> variantAssociatedClinicalTrials();

    @NotNull
    public abstract List<MolecularProfileEvidence> variantLevelEvidences();
}