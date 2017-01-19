package com.hartwig.hmftools.patientreporter.variants;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VariantReport {

    @NotNull
    private final String gene;
    @NotNull
    private final String position;
    @NotNull
    private final String ref;
    @NotNull
    private final String alt;
    @NotNull
    private final String transcript;
    @NotNull
    private final String hgvsCoding;
    @NotNull
    private final String hgvsProtein;
    @NotNull
    private final String consequence;

    @Nullable
    private final String cosmicID;
    private final double alleleFrequency;
    private final int readDepth;

    public VariantReport(@NotNull final String gene, @NotNull final String position, @NotNull final String ref,
            @NotNull final String alt, @NotNull final String transcript, @NotNull final String hgvsCoding,
            @NotNull final String hgvsProtein, @NotNull final String consequence, @Nullable final String cosmicID,
            final double alleleFrequency, final int readDepth) {
        this.gene = gene;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.transcript = transcript;
        this.hgvsCoding = hgvsCoding;
        this.hgvsProtein = hgvsProtein;
        this.consequence = consequence;
        this.cosmicID = cosmicID;
        this.alleleFrequency = alleleFrequency;
        this.readDepth = readDepth;
    }
}
