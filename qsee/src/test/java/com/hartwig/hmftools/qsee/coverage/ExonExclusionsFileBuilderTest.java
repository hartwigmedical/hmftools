package com.hartwig.hmftools.qsee.coverage;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.qsee.common.SampleType;

import org.junit.Test;

public class ExonExclusionsFileBuilderTest
{
    private static final String GENE_NAME = "GENE";
    private static final int SAMPLE_COUNT = 20;
    private static final double MIN_LOW_COVERAGE_SAMPLE_FRACTION = 0.95;

    @Test
    public void canCalcExonExclusionsFromCohort()
    {
        // exon rank, tumor samples with WARN status, normal samples with WARN status, expected excluded
        Object[][] exonRows = {
                { 1, 19, 0, true }, // 19 >= ceil(20 * 0.95)
                { 2, 18, 0, false },
                { 3, 0, 20, false }, // Normal samples ignored
        };

        List<ExonCoverageQc> coverageQcPerExon = new ArrayList<>();

        for(Object[] exonRow : exonRows)
        {
            int exonRank = (int) exonRow[0];
            addCoverageQcPerExon(coverageQcPerExon, SampleType.TUMOR, exonRank, (int) exonRow[1]);
            addCoverageQcPerExon(coverageQcPerExon, SampleType.NORMAL, exonRank, (int) exonRow[2]);
        }

        ExonExclusions exonExclusions = ExonExclusionsFileBuilder.calcExonExclusions(coverageQcPerExon, MIN_LOW_COVERAGE_SAMPLE_FRACTION);

        for(Object[] exonRow : exonRows)
        {
            assertEquals(exonRow[3], exonExclusions.isExcluded(GENE_NAME, (int) exonRow[0]));
        }
    }

    private static void addCoverageQcPerExon(List<ExonCoverageQc> coverageQcPerExon, SampleType sampleType, int exonRank, int warnSampleCount)
    {
        for(int sampleIndex = 0; sampleIndex < SAMPLE_COUNT; ++sampleIndex)
        {
            VariantDetectabilityStatus status = sampleIndex < warnSampleCount ? VariantDetectabilityStatus.WARN : VariantDetectabilityStatus.PASS;

            coverageQcPerExon.add(new ExonCoverageQc(
                    sampleType.name() + "_" + sampleIndex, sampleType,
                    GENE_NAME, exonRank, 100,
                    0.05, status,
                    Map.of()
            ));
        }
    }
}
