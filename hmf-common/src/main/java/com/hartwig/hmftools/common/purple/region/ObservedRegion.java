package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface ObservedRegion extends GenomeRegion {

    int bafCount();

    double observedBAF();

    double observedTumorRatio();

    double observedNormalRatio();
}
