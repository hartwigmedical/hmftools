package com.hartwig.hmftools.common.variant;

class ReadCount {

    private final int alleleReadCount;
    private final int totalReadCount;

    ReadCount(final int alleleReadCount, final int totalReadCount) {
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
