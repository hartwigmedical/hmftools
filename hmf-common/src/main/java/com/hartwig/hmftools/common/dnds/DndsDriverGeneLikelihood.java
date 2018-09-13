package com.hartwig.hmftools.common.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverGeneLikelihood {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract DndsDriverImpactLikelihood missense();

    @NotNull
    public abstract DndsDriverImpactLikelihood nonsense();

    @NotNull
    public abstract DndsDriverImpactLikelihood splice();

    @NotNull
    public abstract DndsDriverImpactLikelihood indel();

}
