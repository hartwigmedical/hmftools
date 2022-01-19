package com.hartwig.hmftools.ckb.json.common;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EvidenceInfo {

    public abstract int id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract MolecularProfileInfo molecularProfile();

    @NotNull
    public abstract TherapyInfo therapy();

    @Nullable
    public abstract IndicationInfo indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<ReferenceInfo> references();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproaches();
}
