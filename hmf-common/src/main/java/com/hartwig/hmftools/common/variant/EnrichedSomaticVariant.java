package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedSomaticVariant implements PurityAdjustedSomaticVariant {

    public abstract String trinucleotideContext();

    public abstract boolean highConfidenceRegion();

    public abstract String microhomology();

    public abstract String repeatSequence();

    public abstract String gene();

    public abstract String effect();

    public abstract int repeatCount();

    public abstract Clonality clonality();
}
