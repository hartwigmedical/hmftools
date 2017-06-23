package com.hartwig.hmftools.common.variant.structural;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantPosition implements GenomePosition {
    public abstract StructuralVariantType type();
}
