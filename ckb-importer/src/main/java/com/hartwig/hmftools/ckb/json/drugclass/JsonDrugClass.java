package com.hartwig.hmftools.ckb.json.drugclass;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonDrugClass implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String drugClass();

    @NotNull
    public abstract LocalDate createDate();

    @NotNull
    public abstract List<DrugInfo> drugs();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproaches();
}
