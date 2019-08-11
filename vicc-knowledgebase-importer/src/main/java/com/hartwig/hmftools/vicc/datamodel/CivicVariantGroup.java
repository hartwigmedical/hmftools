package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariantGroup {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String variants();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String name();
}
