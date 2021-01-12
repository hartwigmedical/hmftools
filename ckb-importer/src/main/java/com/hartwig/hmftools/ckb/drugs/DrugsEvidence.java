package com.hartwig.hmftools.ckb.drugs;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsEvidence {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract DrugsMolecularProfile molecularProfile();

    @NotNull
    public abstract DrugsTherapy therapy();

    @Nullable
    public abstract DrugsIndications indications();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<DrugsEvidenceReferences> references();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();
}
