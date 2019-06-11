package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PhenotypeType {

    @Nullable
    public abstract String source();

    @NotNull
    public abstract String term();

    @NotNull
    public abstract String id();
}
