package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.interpretation.clinicaltrial.ClinicalTrialInterpretation;
import com.hartwig.hmftools.ckb.interpretation.evidence.EvidenceInterpretation;
import com.hartwig.hmftools.ckb.interpretation.knowngenomicalteration.KnownGenomicAlteration;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntryInterpretation {

    public abstract int id();

    public abstract int molecularProfileId();

    @NotNull
    public abstract List<ClinicalTrialInterpretation> clinicalTrialInterpretations();

    @NotNull
    public abstract List<EvidenceInterpretation> evidenceInterpretations();

    @NotNull
    public abstract KnownGenomicAlteration knownAberation();
}
