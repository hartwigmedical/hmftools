package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicDisease {

    @Nullable
    public abstract String doid();

    @NotNull
    public abstract String url();

    @NotNull
    public abstract String displayName();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

}
