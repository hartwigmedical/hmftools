package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public interface Variant extends GenomePosition, AllelicDepth {

    @NotNull
    String gene();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    CodingEffect canonicalCodingEffect();
}
