package com.hartwig.hmftools.ckb.common;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EvidenceInfo {

    @NotNull
    public abstract String id();

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

    @NotNull
    public abstract IndicationInfo indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<ReferenceInfo> reference();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();

    @Nullable
    public abstract List<TreatmentApproach> relevantTreatmentApproach();
}
