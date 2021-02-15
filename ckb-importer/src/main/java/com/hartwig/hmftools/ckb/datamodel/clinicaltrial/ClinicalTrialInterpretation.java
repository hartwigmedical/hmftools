package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.therapyinterpretation.TherapyInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialInterpretation {

    @NotNull
    public abstract ClinicalTrial clinicalTrial();

    @NotNull
    public abstract List<Indication> indications();

    @NotNull
    public abstract List<TherapyInterpretation> therapyInterpretations();


}
