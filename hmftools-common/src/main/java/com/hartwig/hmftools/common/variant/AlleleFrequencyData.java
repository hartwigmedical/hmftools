package com.hartwig.hmftools.common.variant;

class AlleleFrequencyData {

    private final int alleleReadCount;
    private final int totalReadCount;

    AlleleFrequencyData(final int alleleReadCount, final int totalReadCount) {
        this.alleleReadCount = alleleReadCount;
        this.totalReadCount = totalReadCount;
    }

    public int alleleReadCount() {
        return alleleReadCount;
    }

    public int totalReadCount() {
        return totalReadCount;
    }
}
