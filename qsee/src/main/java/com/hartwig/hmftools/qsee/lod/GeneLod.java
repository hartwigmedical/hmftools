package com.hartwig.hmftools.qsee.lod;

import com.hartwig.hmftools.qsee.common.SampleType;

public record GeneLod(
        String sampleId,
        SampleType sampleType,
        String geneName,
        int exonRank,
        double medianDepth,
        double tumorVafLod,
        VariantDetectabilityStatus variantDetectabilityStatus
) {}
