package com.hartwig.hmftools.common.variant;

public interface AllelicDepth {

    int totalReadCount();

    int alleleReadCount();

    default double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }
}
