package com.hartwig.hmftools.qsee.lod;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.qsee.common.SampleType;

import org.junit.Test;

public class GeneLodCalculatorTest
{
    private static final String SAMPLE_ID = "TUMOR_SAMPLE";
    private static final String REPORTABLE_GENE = "REPORTABLE_GENE";
    private static final String NON_REPORTABLE_GENE = "NON_REPORTABLE_GENE";

    private static final Map<DepthRange,Integer> MIXED_DEPTH_COUNTS = Map.of(
            new DepthRange(20, 20), 10,
            new DepthRange(100, 100), 90
    );

    private static final Map<DepthRange,Integer> LOW_DEPTH_COUNTS = Map.of(
            new DepthRange(2, 2), 10,
            new DepthRange(100, 100), 90
    );

    private static final List<ExonDepth> EXON_DEPTHS = List.of(
            new ExonDepth(REPORTABLE_GENE, "1", 100, 199, 1, 100, 92, MIXED_DEPTH_COUNTS),
            new ExonDepth(REPORTABLE_GENE, "1", 300, 399, 2, 100, 90.2, LOW_DEPTH_COUNTS),
            new ExonDepth(NON_REPORTABLE_GENE, "1", 500, 599, 1, 100, 92, MIXED_DEPTH_COUNTS)
    );

    @Test
    public void canCalcGeneLods()
    {
        // purity, exon rank, expected tumor VAF LoD, expected detectable clonal variants
        Object[][] expectedRows = {
                { 0.5, 1, 0.263, VariantDetectabilityStatus.PASS },
                { 0.5, 2, 2.0, VariantDetectabilityStatus.FAIL }, // Depth too low
                { 0.2, 1, 0.658, VariantDetectabilityStatus.WARN_HET_UNRELIABLE },
                { 0.2, 2, 5.0, VariantDetectabilityStatus.FAIL },
                { 0.1, 1, 1.316, VariantDetectabilityStatus.FAIL }, // Purity adjusted VAF > 1
        };

        for(Object[] expectedRow : expectedRows)
        {
            double tumorPurity = (double) expectedRow[0];
            int exonRank = (int) expectedRow[1];

            List<GeneLod> geneLods = GeneLodCalculator.calcGeneLods(
                    SAMPLE_ID, SampleType.TUMOR, EXON_DEPTHS, List.of(REPORTABLE_GENE), tumorPurity);

            assertEquals(2, geneLods.size()); // Non-reportable gene excluded

            GeneLod geneLod = geneLods.stream().filter(x -> x.exonRank() == exonRank).findFirst().orElseThrow();

            assertEquals(REPORTABLE_GENE, geneLod.geneName());
            assertEquals(100, geneLod.medianDepth(), 0);
            assertEquals((double) expectedRow[2], geneLod.tumorVafLod(), 1e-3);
            assertEquals(expectedRow[3], geneLod.variantDetectabilityStatus());
        }
    }
}
