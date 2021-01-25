package com.hartwig.hmftools.ckb.reference;

import java.util.List;

import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReferenceEvidence {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String efficacyEvidence();

    @NotNull
    public abstract ReferenceMolecularProfile molecularProfile();

    @NotNull
    public abstract TherapyInfo therapy();

    @NotNull
    public abstract IndicationInfo indication();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract List<ReferenceReference> reference();

    @NotNull
    public abstract String ampCapAscoEvidenceLevel();

    @NotNull
    public abstract String ampCapAscoInferredTier();

}
