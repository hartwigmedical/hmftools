package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicOrganization {

    @NotNull
    public abstract String url();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract CivicProfileImage profileImage();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String name();
}
