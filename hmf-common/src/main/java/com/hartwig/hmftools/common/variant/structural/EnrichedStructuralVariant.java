package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedStructuralVariant implements StructuralVariant {

    @Nullable
    public abstract Double ploidy();

    @NotNull
    @Override
    public abstract EnrichedStructuralVariantLeg start();

    @NotNull
    @Override
    public abstract EnrichedStructuralVariantLeg end();
}
