package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.purple.segment.PurpleSegmentSource;
import com.hartwig.hmftools.common.region.GenomeRegion;

public interface ObservedRegion extends GenomeRegion {

    PurpleSegmentSource source();

    int bafCount();

    double observedBAF();

    double observedTumorRatio();

    double observedNormalRatio();
}
