package com.hartwig.hmftools.common.variant.structural;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface StructuralVariantLeg extends GenomePosition {

    byte orientation();

    @NotNull
    String homology();

    @Nullable
    Double alleleFrequency();
}
