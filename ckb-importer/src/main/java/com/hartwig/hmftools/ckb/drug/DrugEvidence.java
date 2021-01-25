package com.hartwig.hmftools.ckb.drug;

import java.util.List;

import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugEvidence {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract DrugMolecularProfile molecularProfile();

    @NotNull
    public abstract TherapyInfo therapy();

    @Nullable
    public abstract DrugIndication indications();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<DrugEvidenceReference> references();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();
}
