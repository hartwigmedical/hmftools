package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;

public interface ObservedRegion extends GenomeRegion {

    boolean ratioSupport();

    StructuralVariantSupport structuralVariantSupport();

    int bafCount();

    double observedBAF();

    int observedTumorRatioCount();

    double observedTumorRatio();

    double observedNormalRatio();

    FreecStatus status();
}
