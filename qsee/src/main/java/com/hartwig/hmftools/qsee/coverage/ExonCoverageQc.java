package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.feature.FeatureType.EXON_VARIANT_LOD;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.status.QcStatusType;
import com.hartwig.hmftools.qsee.status.ThresholdRegistry;

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
    public static final String LOD_VAF = "LOD_VAF";
    public static final String MEDIAN_DEPTH = "MEDIAN_DEPTH";

    private static VariantDetectabilityStatus determineVariantDetectabilityStatus(
            SampleType sampleType,
            String geneName,
            int exonRank,
            double medianDepth,
            double lodVaf,
            ExonExclusions exonExclusions,
            ThresholdRegistry qcThresholds
    ){
        if(exonExclusions.isExcluded(geneName, exonRank))
            return VariantDetectabilityStatus.EXCLUDED;

        boolean lodVafFails = isWarn(qcThresholds, sampleType, LOD_VAF, lodVaf);
        boolean medianDepthFails = isWarn(qcThresholds, sampleType, MEDIAN_DEPTH, medianDepth);

        if(lodVafFails || medianDepthFails)
            return VariantDetectabilityStatus.WARN;

        return VariantDetectabilityStatus.PASS;
    }

    private static boolean isWarn(ThresholdRegistry qcThresholds, SampleType sampleType, String featureName, double value)
    {
        return qcThresholds.getThreshold(sampleType, EXON_VARIANT_LOD, featureName, QcStatusType.WARN)
                .getSampleQcStatus(value).type() == QcStatusType.WARN;
    }

    public static ExonCoverageQc from(
            String sampleId,
            SampleType sampleType,
            ExonDepth exonDepth,
            ExonExclusions exonExclusions,
            ThresholdRegistry qcThresholds,
            int lodMinAltReads
    ){
        double lodVaf = LimitOfDetectionCalcs.calcMinSampleVafToDetectVariant(exonDepth.DepthCounts, lodMinAltReads);

        VariantDetectabilityStatus variantDetectabilityStatus = determineVariantDetectabilityStatus(
                sampleType,
                exonDepth.Gene, exonDepth.ExonRank, exonDepth.MedianDepth,
                lodVaf, exonExclusions, qcThresholds
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
            List<String> reportableGenes,
            ThresholdRegistry qcThresholds,
            int lodMinAltReads
    ){
        List<ExonCoverageQc> coverageQcPerExon = new ArrayList<>();

        for(ExonDepth exonDepth : exonDepths)
        {
            if(!reportableGenes.contains(exonDepth.Gene))
                continue;

            coverageQcPerExon.add(from(sampleId, sampleType, exonDepth, exonExclusions, qcThresholds, lodMinAltReads));
        }

        return coverageQcPerExon;
    }

    public static double calcLowCoverageExonProportion(List<ExonCoverageQc> coverageQcPerExon)
    {
        long warnCount = coverageQcPerExon.stream().filter(x -> x.variantDetectabilityStatus() == VariantDetectabilityStatus.WARN).count();
        long passCount = coverageQcPerExon.stream().filter(x -> x.variantDetectabilityStatus() == VariantDetectabilityStatus.PASS).count();

        long nonExcludedCount = warnCount + passCount;
        return nonExcludedCount > 0 ? (double) warnCount / nonExcludedCount : 0;
    }
}
