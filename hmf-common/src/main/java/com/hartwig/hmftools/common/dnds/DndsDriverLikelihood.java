package com.hartwig.hmftools.common.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DndsDriverLikelihood {

    public abstract String gene();

    public abstract double missenseUnadjustedDriverLikelihood();

    public abstract double missenseProbabilityDriver();

    public abstract double missenseProbabilityVariantNonDriverFactor();
    
}
