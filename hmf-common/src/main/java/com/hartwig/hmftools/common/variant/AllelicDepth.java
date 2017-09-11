package com.hartwig.hmftools.common.variant;

public interface AllelicDepth {

    int totalReadCount();

    int alleleReadCount();

    double alleleFrequency();
}
