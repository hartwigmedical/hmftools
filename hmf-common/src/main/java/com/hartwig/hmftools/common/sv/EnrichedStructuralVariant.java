package com.hartwig.hmftools.common.sv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedStructuralVariant implements StructuralVariant {

    @Nullable
    public abstract Double junctionCopyNumber();

    @NotNull
    @Override
    public abstract EnrichedStructuralVariantLeg start();

    @Nullable
    @Override
    public abstract EnrichedStructuralVariantLeg end();
}
