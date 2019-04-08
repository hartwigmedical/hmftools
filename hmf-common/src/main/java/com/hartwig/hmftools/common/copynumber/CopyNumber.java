package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface CopyNumber extends GenomeRegion {

    int value();

    default CopyNumberAlteration alteration() {
        return CopyNumberAlteration.fromCopyNumber(value());
    }
}
