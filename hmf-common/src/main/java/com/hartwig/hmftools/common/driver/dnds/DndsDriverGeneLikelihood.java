package com.hartwig.hmftools.common.driver.dnds;

import com.hartwig.hmftools.common.driver.DriverImpact;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverGeneLikelihood
{
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

    @NotNull
    public DndsDriverImpactLikelihood select(final DriverImpact impact)
    {
        switch(impact)
        {
            case NONSENSE:
                return nonsense();
            case MISSENSE:
                return missense();
            case INFRAME:
            case FRAMESHIFT:
                return indel();
            case SPLICE:
                return splice();
            default:
                throw new IllegalArgumentException();
        }
    }
}
