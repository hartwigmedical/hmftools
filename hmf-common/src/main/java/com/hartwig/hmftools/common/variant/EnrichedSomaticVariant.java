package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedSomaticVariant implements PurityAdjustedSomaticVariant {

    @NotNull
    public abstract String trinucleotideContext();

    public abstract boolean highConfidenceRegion();

    @NotNull
    public abstract String microhomology();

    @NotNull
    public abstract String repeatSequence();

    public abstract int repeatCount();

    @NotNull
    public abstract Clonality clonality();
}
