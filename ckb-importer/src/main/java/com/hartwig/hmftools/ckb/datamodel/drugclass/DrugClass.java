package com.hartwig.hmftools.ckb.datamodel.drugclass;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbJsonObject;
import com.hartwig.hmftools.ckb.datamodel.common.DrugInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClass implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String drugClass();

    @Nullable
    public abstract Date createDate();

    @NotNull
    public abstract List<DrugInfo> drug();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproach();
}
