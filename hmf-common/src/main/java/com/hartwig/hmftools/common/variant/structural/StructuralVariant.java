package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariant {

    public abstract String id();

    @Nullable
    public abstract String mateId();

    public abstract String startChromosome();

    public abstract String endChromosome();

    public abstract long startPosition();

    public abstract long endPosition();

    public abstract byte startOrientation();

    public abstract byte endOrientation();

    public abstract String startHomology();

    public abstract String endHomology();

    @Nullable
    public abstract Double startAF();

    @Nullable
    public abstract Double endAF();

    public abstract String insertSequence();

    public abstract StructuralVariantType type();
}
