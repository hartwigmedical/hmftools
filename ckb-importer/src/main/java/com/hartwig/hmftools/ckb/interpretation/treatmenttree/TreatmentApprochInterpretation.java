package com.hartwig.hmftools.ckb.interpretation.treatmenttree;

import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApprochInterpretation {

    @Nullable
    public abstract DrugClass drugClass();

    @Nullable
    public abstract Therapy therapy();
}
