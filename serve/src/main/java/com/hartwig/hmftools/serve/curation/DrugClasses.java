package com.hartwig.hmftools.serve.curation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClasses {

    @NotNull
    public abstract String drugClass();

    @NotNull
    public abstract String curatedDrugClass();
}