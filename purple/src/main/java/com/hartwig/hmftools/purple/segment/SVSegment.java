package com.hartwig.hmftools.purple.segment;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class SVSegment implements GenomePosition {

    @NotNull
    public abstract StructuralVariantType type();
}
