package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface CopyNumber extends GenomeRegion {

    int NORMAL_HUMAN_COPY_NUMBER = 2;

    int value();

    default boolean isGain() {
        return value() > NORMAL_HUMAN_COPY_NUMBER;
    }

    default boolean isLoss() {
        return value() < NORMAL_HUMAN_COPY_NUMBER;
    }
}
