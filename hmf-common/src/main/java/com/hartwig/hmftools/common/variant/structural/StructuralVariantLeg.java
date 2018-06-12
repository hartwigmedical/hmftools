package com.hartwig.hmftools.common.variant.structural;

import com.hartwig.hmftools.common.position.GenomeInterval;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface StructuralVariantLeg extends GenomeInterval {

    byte orientation();

    @NotNull
    String homology();

    @Nullable
    Double alleleFrequency();

    @NotNull
    long impreciseHomologyIntervalStart();

    @NotNull
    long impreciseHomologyIntervalEnd();
}
