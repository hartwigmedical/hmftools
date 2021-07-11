package com.hartwig.hmftools.common.sv;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface EnrichedStructuralVariantLeg extends StructuralVariantLeg {

    @Nullable
    Double adjustedAlleleFrequency();

    @Nullable
    Double adjustedCopyNumber();

    @Nullable
    Double adjustedCopyNumberChange();

    @Nullable
    String refGenomeContext();
}
