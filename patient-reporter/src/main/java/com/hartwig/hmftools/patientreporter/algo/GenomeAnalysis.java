package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.purple.PurpleAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;

import org.jetbrains.annotations.NotNull;

public class GenomeAnalysis {
    @NotNull
    private final String sample;
    @NotNull
    private final VariantAnalysis variantAnalysis;
    @NotNull
    private final CopyNumberAnalysis copyNumberAnalysis;
    @NotNull
    private final PurpleAnalysis purpleAnalysis;

    GenomeAnalysis(@NotNull final String sample, @NotNull final VariantAnalysis variantAnalysis,
            @NotNull final CopyNumberAnalysis copyNumberAnalysis, @NotNull final PurpleAnalysis purpleAnalysis) {
        this.sample = sample;
        this.variantAnalysis = variantAnalysis;
        this.copyNumberAnalysis = copyNumberAnalysis;
        this.purpleAnalysis = purpleAnalysis;
    }

    @NotNull
    PurpleAnalysis purpleAnalysis() {
        return purpleAnalysis;
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
