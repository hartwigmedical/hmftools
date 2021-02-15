package com.hartwig.hmftools.ckb.datamodel;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialInterpretation;
import com.hartwig.hmftools.ckb.datamodel.evidence.EvidenceInterpretation;
import com.hartwig.hmftools.ckb.datamodel.knowngenomicalteration.KnownGenomicAlteration;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntry {

    public abstract int profileId();

    @NotNull
    public abstract String profileName();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract List<ClinicalTrialInterpretation> clinicalTrialInterpretations();

    @NotNull
    public abstract List<EvidenceInterpretation> evidenceInterpretations();

    @NotNull
    public abstract KnownGenomicAlteration knownGenomicAlteration();
}