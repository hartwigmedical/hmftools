package com.hartwig.hmftools.patientdb.data;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticVariantData {
    @NotNull
    private final String gene;
    @NotNull
    private final String position;
    @NotNull
    private final String ref;
    @NotNull
    private final String alt;
    @Nullable
    private final String cosmicID;
    private final int totalReadCount;
    private final int alleleReadCount;

    public SomaticVariantData(@NotNull final String gene, @NotNull final String position, @NotNull final String ref,
            @NotNull final String alt, @Nullable final String cosmicID, final int totalReadCount,
            final int alleleReadCount) {
        this.gene = gene;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.cosmicID = cosmicID;
        this.totalReadCount = totalReadCount;
        this.alleleReadCount = alleleReadCount;
    }

    @NotNull
    public String gene() {
        return gene;
    }

    @NotNull
    public String position() {
        return position;
    }

    @NotNull
    public String ref() {
        return ref;
    }

    @NotNull
    public String alt() {
        return alt;
    }

    @Nullable
    public String cosmicID() {
        return cosmicID;
    }

    public int alleleReadCount() {
        return alleleReadCount;
    }

    public int totalReadCount() {
        return totalReadCount;
    }
}
