package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariantTypes {

    @NotNull
    public abstract String displayName();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String url();

    @NotNull
    public abstract String soId();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();
}
