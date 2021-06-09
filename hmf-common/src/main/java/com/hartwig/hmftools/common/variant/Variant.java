package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public interface Variant extends GenomePosition, AllelicDepth {

    @NotNull
    VariantType type();

    @NotNull
    String gene();

    @NotNull
    String ref();

    @NotNull
    String alt();

    @NotNull
    String canonicalTranscript();

    @NotNull
    CodingEffect canonicalCodingEffect();

    @NotNull
    String canonicalHgvsCodingImpact();

    @NotNull
    String canonicalHgvsProteinImpact();

    @NotNull
    default String genomicEvent() {
        String description = canonicalCodingEffect() == CodingEffect.SPLICE ? canonicalHgvsCodingImpact() : canonicalHgvsProteinImpact();
        return this.gene() + " " + description;
    }
}
