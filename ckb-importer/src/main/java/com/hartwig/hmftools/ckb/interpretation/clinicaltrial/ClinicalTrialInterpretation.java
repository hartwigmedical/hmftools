package com.hartwig.hmftools.ckb.interpretation.clinicaltrial;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialInterpretation {

    @NotNull
    public abstract ClinicalTrial clinicalTrials();

    @NotNull
    public abstract List<Indication> indications();

    @NotNull
    public abstract List<TherapyInterpretation> therapyInterpretations();


}
