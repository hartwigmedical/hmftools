package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.feature.FeatureType.EXON_VARIANT_LOD;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.status.QcStatusType;
import com.hartwig.hmftools.qsee.status.QcThreshold;
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

    private static VariantDetectabilityStatus determineVariantDetectabilityStatus(
            SampleType sampleType,
            String geneName,
            int exonRank,
            double lodVaf,
            ExonExclusions exonExclusions,
            ThresholdRegistry qcThresholds
    ){
        if(exonExclusions.isExcluded(geneName, exonRank))
            return VariantDetectabilityStatus.EXCLUDED;

        QcThreshold threshold = qcThresholds.getThreshold(sampleType, EXON_VARIANT_LOD, LOD_VAF, QcStatusType.WARN);
        QcStatusType lodVafStatus = threshold.getSampleQcStatus(lodVaf).type();

        if(lodVafStatus == QcStatusType.WARN)
            return VariantDetectabilityStatus.WARN;

        return VariantDetectabilityStatus.PASS;
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
                exonDepth.Gene, exonDepth.ExonRank,
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
