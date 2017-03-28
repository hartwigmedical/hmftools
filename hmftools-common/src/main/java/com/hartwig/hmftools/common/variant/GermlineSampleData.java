package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

public class GermlineSampleData {

    @NotNull
    private final String originalSampleData;
    @NotNull
    private final String genoType;
    private final int totalReadCount;
    private final int alleleReadCount;

    GermlineSampleData(@NotNull final String originalSampleData, @NotNull final String genoType,
            final int totalReadCount, final int alleleReadCount) {
        this.originalSampleData = originalSampleData;
        this.genoType = genoType;
        this.totalReadCount = totalReadCount;
        this.alleleReadCount = alleleReadCount;
    }
}
