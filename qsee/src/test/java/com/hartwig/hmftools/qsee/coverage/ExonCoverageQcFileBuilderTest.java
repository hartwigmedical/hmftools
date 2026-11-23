package com.hartwig.hmftools.qsee.coverage;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.qsee.common.SampleType;

import org.junit.Test;

public class ExonCoverageQcFileBuilderTest
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

    private static final Map<DepthRange,Integer> HIGH_DEPTH_COUNTS = Map.of(
            new DepthRange(500, 500), 100
    );

    private static final List<ExonDepth> EXON_DEPTHS = List.of(
            new ExonDepth(REPORTABLE_GENE, "1", 100, 199, 1, 100, 92, MIXED_DEPTH_COUNTS),
            new ExonDepth(REPORTABLE_GENE, "1", 300, 399, 2, 100, 90.2, LOW_DEPTH_COUNTS),
            new ExonDepth(REPORTABLE_GENE, "1", 500, 599, 3, 100, 500, HIGH_DEPTH_COUNTS),
            new ExonDepth(NON_REPORTABLE_GENE, "1", 500, 599, 1, 100, 92, MIXED_DEPTH_COUNTS)
    );

    @Test
    public void canCalcCoverageQcPerExon()
    {
        // exon rank, expected LoD VAF, expected status
        Object[][] expectedRows = {
                { 1, 0.132, VariantDetectabilityStatus.WARN },
                { 2, 1.0, VariantDetectabilityStatus.WARN }, // Depth too low
                { 3, 0.0125, VariantDetectabilityStatus.PASS },
        };

        List<ExonCoverageQc> coverageQcPerExon = ExonCoverageQcFileBuilder.calcCoverageQcPerExon(SAMPLE_ID, SampleType.TUMOR, EXON_DEPTHS, List.of(REPORTABLE_GENE));

        assertEquals(3, coverageQcPerExon.size()); // Non-reportable gene excluded

        for(Object[] expectedRow : expectedRows)
        {
            int exonRank = (int) expectedRow[0];

            ExonCoverageQc exonCoverageQc = coverageQcPerExon.stream().filter(x -> x.exonRank() == exonRank).findFirst().orElseThrow();

            assertEquals(REPORTABLE_GENE, exonCoverageQc.geneName());
            assertEquals(100, exonCoverageQc.medianDepth(), 0);
            assertEquals((double) expectedRow[1], exonCoverageQc.lodVaf(), 1e-3);
            assertEquals(expectedRow[2], exonCoverageQc.variantDetectabilityStatus());
        }
    }
}
