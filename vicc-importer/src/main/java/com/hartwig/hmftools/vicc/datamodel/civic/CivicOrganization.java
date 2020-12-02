package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicOrganization {

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String url();

    @Nullable
    public abstract CivicProfileImage profileImage();

    @Nullable
    public abstract String id();

    @Nullable
    public abstract String description();
}
