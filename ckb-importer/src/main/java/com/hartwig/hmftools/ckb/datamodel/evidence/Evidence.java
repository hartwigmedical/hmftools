package com.hartwig.hmftools.ckb.datamodel.evidence;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproaches.RelevantTreatmentApproaches;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Evidence {

    public abstract int id();

    @Nullable
    public abstract Therapy therapy();

    @NotNull
    public abstract Indication indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<RelevantTreatmentApproaches> relevantTreatmentApproaches();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();

    @NotNull
    public abstract List<Reference> references();
}
