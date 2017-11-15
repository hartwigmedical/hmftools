package com.hartwig.hmftools.common.purple.copynumber.sv;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class StructuralVariantLeg implements GenomePosition {

    public abstract int orientation();

    public abstract double vaf();
}
