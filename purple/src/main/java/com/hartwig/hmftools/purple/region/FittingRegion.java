package com.hartwig.hmftools.purple.region;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;

public interface FittingRegion extends GenomeRegion
{
    GermlineStatus germlineStatus();
    int bafCount();
    double observedBAF();
    double observedNormalRatio();
    double observedTumorRatio();
}
