package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface CopyNumber extends GenomeRegion {

    int value();

    default CopyNumberAlteration alteration() {
        return CopyNumberAlteration.fromCopyNumber(value());
    }

    default boolean isGain() {
        return alteration().equals(CopyNumberAlteration.GAIN);
    }

    default boolean isLoss() {
        return alteration().equals(CopyNumberAlteration.LOSS);
    }
}
