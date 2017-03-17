package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;

import org.jetbrains.annotations.NotNull;

public class GenomeAnalysis {
    @NotNull
    private final String sample;
    @NotNull
    private final VariantAnalysis variantAnalysis;
    @NotNull
    private final CopyNumberAnalysis copyNumberAnalysis;

    GenomeAnalysis(@NotNull final String sample, @NotNull final VariantAnalysis variantAnalysis,
            @NotNull final CopyNumberAnalysis copyNumberAnalysis) {
        this.sample = sample;
        this.variantAnalysis = variantAnalysis;
        this.copyNumberAnalysis = copyNumberAnalysis;
    }

    @NotNull
    public String sample() {
        return sample;
    }

    @NotNull
    public VariantAnalysis variantAnalysis() {
        return variantAnalysis;
    }

    @NotNull
    public CopyNumberAnalysis copyNumberAnalysis() {
        return copyNumberAnalysis;
    }
}
