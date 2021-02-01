package com.hartwig.hmftools.ckb.datamodel.drugclass;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.DrugInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClass {

    public abstract int id();

    @NotNull
    public abstract String drugClass();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<DrugInfo> drug();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproach();
}
