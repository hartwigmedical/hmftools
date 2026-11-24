package com.hartwig.hmftools.qsee.coverage;

import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.qsee.common.SampleType;

public record ExonCoverageQc(
        String sampleId,
        SampleType sampleType,
        String geneName,
        int exonRank,
        double medianDepth,
        double sampleVafLod,
        double tumorVafLod,
        VariantDetectabilityStatus variantDetectabilityStatus,
        Map<DepthRange,Integer> depthCounts
) {}
