package com.hartwig.hmftools.ckb.drugClasses;

import java.util.List;

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
    public abstract List<DrugClassDrugs> drugs();

    @NotNull
    public abstract List<DrugClassTreatmentApproaches> treatmentApproaches();
}
