package com.hartwig.hmftools.common.variant.structural;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;

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
    CopyNumberMethod method();

    @Nullable
    String refGenomeContext();
}
