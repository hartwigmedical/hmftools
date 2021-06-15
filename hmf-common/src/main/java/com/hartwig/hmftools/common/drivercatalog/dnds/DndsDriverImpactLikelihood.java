package com.hartwig.hmftools.common.drivercatalog.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverImpactLikelihood
{
    public abstract double driversPerSample();

    public abstract double passengersPerMutation();

}
