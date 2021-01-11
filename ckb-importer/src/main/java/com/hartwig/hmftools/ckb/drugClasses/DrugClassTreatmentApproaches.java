package com.hartwig.hmftools.ckb.drugClasses;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClassTreatmentApproaches {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String profileName();
}
