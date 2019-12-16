package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicDrug {

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String pubchemId();

    @NotNull
    public abstract String id();
}
