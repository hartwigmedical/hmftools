package com.hartwig.hmftools.ckb.interpretation.clinicaltrialtree;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;

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


}
