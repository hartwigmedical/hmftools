package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.coverage.LimitOfDetectionCalcs.LOD_VAF_WARN_THRESHOLD;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.qsee.common.SampleType;

public record ExonCoverageQc(
        String sampleId,
        SampleType sampleType,
        String geneName,
        int exonRank,
        double medianDepth,
        double lodVaf,
        VariantDetectabilityStatus variantDetectabilityStatus,
        Map<DepthRange,Integer> depthCounts
)
{
    private static VariantDetectabilityStatus determineVariantDetectabilityStatus(
            double lodVaf,
            String geneName,
            int exonRank,
            ExonExclusions exonExclusions
    ){
        if(exonExclusions.isExcluded(geneName, exonRank))
            return VariantDetectabilityStatus.EXCLUDED;

        if(lodVaf > LOD_VAF_WARN_THRESHOLD)
            return VariantDetectabilityStatus.WARN;

        return VariantDetectabilityStatus.PASS;
    }

    public static ExonCoverageQc from(String sampleId, SampleType sampleType, ExonDepth exonDepth, ExonExclusions exonExclusions)
    {
        double lodVaf = LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant(exonDepth.DepthCounts);

        VariantDetectabilityStatus variantDetectabilityStatus = determineVariantDetectabilityStatus(
                lodVaf,
                exonDepth.Gene, exonDepth.ExonRank, exonExclusions
        );

        return new ExonCoverageQc(
                sampleId, sampleType,
                exonDepth.Gene, exonDepth.ExonRank, exonDepth.MedianDepth,
                lodVaf, variantDetectabilityStatus,
                exonDepth.DepthCounts
        );
    }

    public static List<ExonCoverageQc> fromExonDepths(
            String sampleId,
            SampleType sampleType,
            List<ExonDepth> exonDepths,
            ExonExclusions exonExclusions,
            List<String> reportableGenes
    ){
        List<ExonCoverageQc> coverageQcPerExon = new ArrayList<>();

        for(ExonDepth exonDepth : exonDepths)
        {
            if(!reportableGenes.contains(exonDepth.Gene))
                continue;

            coverageQcPerExon.add(from(sampleId, sampleType, exonDepth, exonExclusions));
        }

        return coverageQcPerExon;
    }
}
