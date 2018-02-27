package com.hartwig.hmftools.common.variant;

class AllelicDepthImpl implements AllelicDepth {

    private final int alleleReadCount;
    private final int totalReadCount;

    AllelicDepthImpl(final int alleleReadCount, final int totalReadCount) {
        this.alleleReadCount = alleleReadCount;
        this.totalReadCount = totalReadCount;
    }

    @Override
    public int totalReadCount() {
        return totalReadCount;
    }

    @Override
    public int alleleReadCount() {
        return alleleReadCount;
    }
}
