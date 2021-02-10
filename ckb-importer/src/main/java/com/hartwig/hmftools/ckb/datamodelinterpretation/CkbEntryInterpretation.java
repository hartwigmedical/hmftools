package com.hartwig.hmftools.ckb.datamodelinterpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.interpretation.clinicaltrial.ClinicalTrialInterpretation;
import com.hartwig.hmftools.ckb.interpretation.evidence.EvidenceInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CkbEntryInterpretation {

    public abstract int id();

    @NotNull
    public abstract List<ClinicalTrialInterpretation> clinicalTrialInterpretation();

    @NotNull
    public abstract List<EvidenceInterpretation> evidenceInterpretation();
}
