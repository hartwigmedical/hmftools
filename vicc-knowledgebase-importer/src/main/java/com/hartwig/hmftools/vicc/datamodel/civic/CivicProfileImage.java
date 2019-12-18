package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicProfileImage {

    @NotNull
    public abstract String x14();

    @NotNull
    public abstract String x32();

    @NotNull
    public abstract String x64();

    @NotNull
    public abstract String x128();

    @NotNull
    public abstract String x256();
}
