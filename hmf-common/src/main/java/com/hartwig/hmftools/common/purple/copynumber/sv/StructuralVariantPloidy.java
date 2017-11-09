package com.hartwig.hmftools.common.purple.copynumber.sv;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class StructuralVariantPloidy implements GenomePosition {

    public abstract double unweightedImpliedPloidy();

    public abstract double averageImpliedPloidy();

    public abstract double weight();

    public abstract int orientation();

    public abstract boolean alternate();

    public abstract double adjacentCopyNumber();
}
