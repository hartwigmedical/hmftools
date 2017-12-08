package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;

import org.jetbrains.annotations.NotNull;

class GenomeAnalysis {
    @NotNull
    private final String sample;
    @NotNull
    private final VariantAnalysis variantAnalysis;
    @NotNull
    private final PurpleAnalysis purpleAnalysis;
    @NotNull
    private final StructuralVariantAnalysis structuralVariantAnalysis;

    GenomeAnalysis(@NotNull final String sample, @NotNull final VariantAnalysis variantAnalysis,
            @NotNull final PurpleAnalysis purpleAnalysis, @NotNull final StructuralVariantAnalysis structuralVariantAnalysis) {
        this.sample = sample;
        this.variantAnalysis = variantAnalysis;
        this.purpleAnalysis = purpleAnalysis;
        this.structuralVariantAnalysis = structuralVariantAnalysis;
    }

    @NotNull
    PurpleAnalysis purpleAnalysis() {
        return purpleAnalysis;
    }

    @NotNull
    String sample() {
        return sample;
    }

    @NotNull
    VariantAnalysis variantAnalysis() {
        return variantAnalysis;
    }

    @NotNull
    StructuralVariantAnalysis structuralVariantAnalysis() {
        return structuralVariantAnalysis;
    }
}
