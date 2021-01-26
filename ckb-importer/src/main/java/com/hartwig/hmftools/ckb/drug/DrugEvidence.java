package com.hartwig.hmftools.ckb.drug;

import java.util.List;

import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.ReferenceInfo;
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
    public abstract MolecularProfileInfo molecularProfile();

    @NotNull
    public abstract TherapyInfo therapy();

    @Nullable
    public abstract IndicationInfo indications();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<ReferenceInfo> references();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();
}
