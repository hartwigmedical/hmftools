package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariant {

    public abstract String startChromosome();

    public abstract String endChromosome();

    public abstract long startPosition();

    public abstract long endPosition();

    public abstract StructuralVariantType type();
}
