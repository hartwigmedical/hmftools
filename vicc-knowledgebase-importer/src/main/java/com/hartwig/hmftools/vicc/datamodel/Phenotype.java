package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Phenotype {

    @NotNull
    public abstract PhenotypeType type();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String family();

    @NotNull
    public abstract String id();
}
