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
    private static final double MIN_HIGH_LOD_VAF_SAMPLE_FRACTION = 0.95;

    private static final double HIGH_LOD_VAF = 0.10;
    private static final double LOW_LOD_VAF = 0.02;

    @Test
    public void canCalcExonExclusionsFromCohort()
    {
        // exon rank, tumor samples with high LoD VAF, normal samples with high LoD VAF, expected excluded
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

        ExonExclusions exonExclusions = ExonExclusionsFileBuilder.calcExonExclusions(coverageQcPerExon, MIN_HIGH_LOD_VAF_SAMPLE_FRACTION);

        for(Object[] exonRow : exonRows)
        {
            assertEquals(exonRow[3], exonExclusions.isExcluded(GENE_NAME, (int) exonRow[0]));
        }
    }

    private static void addCoverageQcPerExon(List<ExonCoverageQc> coverageQcPerExon, SampleType sampleType, int exonRank, int highLodVafSampleCount)
    {
        for(int sampleIndex = 0; sampleIndex < SAMPLE_COUNT; ++sampleIndex)
        {
            double lodVaf = sampleIndex < highLodVafSampleCount ? HIGH_LOD_VAF : LOW_LOD_VAF;

            coverageQcPerExon.add(new ExonCoverageQc(
                    sampleType.name() + "_" + sampleIndex, sampleType,
                    GENE_NAME, exonRank, 100,
                    lodVaf, VariantDetectabilityStatus.PASS,
                    Map.of()
            ));
        }
    }
}
