package com.hartwig.hmftools.ckb.drugclass;

import java.util.List;

import com.hartwig.hmftools.ckb.common.DrugInfo;
import com.hartwig.hmftools.ckb.common.TreatmentApproachInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClass {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String drugClass();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<DrugInfo> drug();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproach();
}
