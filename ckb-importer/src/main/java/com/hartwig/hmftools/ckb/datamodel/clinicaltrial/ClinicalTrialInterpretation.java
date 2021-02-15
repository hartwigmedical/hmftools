package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.Indication;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.Therapy;

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
    public abstract List<Therapy> therapies();


}
