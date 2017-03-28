package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public class GermlineSampleData {

    @NotNull
    private final String genoType;
    private final int totalReadCount;
    private final int alleleReadCount;

    public GermlineSampleData(@NotNull final String genoType, final int totalReadCount, final int alleleReadCount) {
        this.genoType = genoType;
        this.totalReadCount = totalReadCount;
        this.alleleReadCount = alleleReadCount;
    }

    @NotNull
    public String genoType() {
        return genoType;
    }

    public int totalReadCount() {
        return totalReadCount;
    }

    public int alleleReadCount() {
        return alleleReadCount;
    }

    public double alleleFrequency() {
        return (double) alleleReadCount / totalReadCount;
    }
}
