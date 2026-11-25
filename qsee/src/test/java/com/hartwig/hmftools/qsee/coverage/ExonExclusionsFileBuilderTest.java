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
    private static final double MIN_LOW_COVERAGE_LOD_VAF = 0.07;
    private static final double MIN_LOW_COVERAGE_COHORT_FRACTION = 0.3;

    private static final double LOW_LOD_VAF = 0.01;

    @Test
    public void canCalcExonExclusionsFromCohort()
    {
        // exon rank, tumor samples with high LoD VAF, high LoD VAF, normal samples with high LoD VAF, expected excluded
        Object[][] exonRows = {
                { 1, 6, 0.07, 0, true }, // 6 >= ceil(20 * 0.3), 0.07 >= 0.07
                { 2, 5, 0.07, 0, false },
                { 3, 20, 0.069, 0, false },
                { 4, 0, 0.07, 20, false }, // Normal samples ignored
        };

        List<ExonCoverageQc> coverageQcPerExon = new ArrayList<>();

        for(Object[] exonRow : exonRows)
        {
            int exonRank = (int) exonRow[0];
            double highLodVaf = (double) exonRow[2];
            addCoverageQcPerExon(coverageQcPerExon, SampleType.TUMOR, exonRank, (int) exonRow[1], highLodVaf);
            addCoverageQcPerExon(coverageQcPerExon, SampleType.NORMAL, exonRank, (int) exonRow[3], highLodVaf);
        }

        ExonExclusions exonExclusions = ExonExclusionsFileBuilder.calcExonExclusions(
                coverageQcPerExon, MIN_LOW_COVERAGE_LOD_VAF, MIN_LOW_COVERAGE_COHORT_FRACTION);

        for(Object[] exonRow : exonRows)
        {
            assertEquals(exonRow[4], exonExclusions.isExcluded(GENE_NAME, (int) exonRow[0]));
        }
    }

    private static void addCoverageQcPerExon(
            List<ExonCoverageQc> coverageQcPerExon, SampleType sampleType, int exonRank, int highLodSampleCount, double highLodVaf)
    {
        for(int sampleIndex = 0; sampleIndex < SAMPLE_COUNT; ++sampleIndex)
        {
            double lodVaf = sampleIndex < highLodSampleCount ? highLodVaf : LOW_LOD_VAF;

            coverageQcPerExon.add(new ExonCoverageQc(
                    sampleType.name() + "_" + sampleIndex, sampleType,
                    GENE_NAME, exonRank, 100,
                    lodVaf, VariantDetectabilityStatus.PASS,
                    Map.of()
            ));
        }
    }
}
