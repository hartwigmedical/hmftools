package com.hartwig.hmftools.ckb.interpretation.evidence;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ReferenceExtend;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication;
import com.hartwig.hmftools.ckb.interpretation.common.molecularprofileinterpretation.MolecularProfileInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EvidenceInterpretation {

    public abstract int id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract MolecularProfileInterpretation variantInterpretation();

    @Nullable
    public abstract TherapyInterpretation therapyInterpretation();

    @NotNull
    public abstract Indication indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<ReferenceExtend> references();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();
}
