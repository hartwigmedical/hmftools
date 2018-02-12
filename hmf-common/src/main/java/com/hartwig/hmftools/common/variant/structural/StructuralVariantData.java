package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
public abstract class StructuralVariantData {

    @Nullable
    public abstract String id();

    public abstract String startChromosome();
    public abstract String endChromosome();
    public abstract long startPosition();
    public abstract long endPosition();

    public abstract byte startOrientation();
    public abstract byte endOrientation();
    public abstract StructuralVariantType type();
}
