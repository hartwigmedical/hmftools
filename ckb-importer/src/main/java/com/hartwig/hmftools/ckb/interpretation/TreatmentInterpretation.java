package com.hartwig.hmftools.ckb.interpretation;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentInterpretation {

    @NotNull
    public abstract TreatmentApproach treatmentApproach();

    @NotNull
    public abstract DrugClass drugClass();

    @NotNull
    public abstract Therapy therapy();

    @NotNull
    public abstract Drug drug();
}
