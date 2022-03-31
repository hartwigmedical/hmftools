package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface Variant extends GenomePosition, AllelicDepth {

    @NotNull
    VariantType type();

    @NotNull
    String gene();

    @Nullable
    String transcript();

    @Nullable
    Boolean isCanonical();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    String canonicalTranscript();

    @NotNull
    String canonicalEffect();

    @NotNull
    CodingEffect canonicalCodingEffect();

    @NotNull
    String canonicalHgvsCodingImpact();

    @NotNull
    String canonicalHgvsProteinImpact();
}
