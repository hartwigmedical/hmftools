package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;

public interface ObservedRegion extends GenomeRegion {

    boolean ratioSupport();

    SegmentSupport support();

    int bafCount();

    double observedBAF();

    int observedTumorRatioCount();

    double observedTumorRatio();

    double observedNormalRatio();

    GermlineStatus status();

    boolean svCluster();

    double gcContent();
}
