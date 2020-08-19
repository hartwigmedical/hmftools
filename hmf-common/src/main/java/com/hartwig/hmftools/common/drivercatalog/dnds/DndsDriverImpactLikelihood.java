package com.hartwig.hmftools.common.drivercatalog.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverImpactLikelihood {

    @Deprecated
    public abstract double dndsLikelihood();

    public abstract double pDriver();

    public abstract double pVariantNonDriverFactor();
    
}
