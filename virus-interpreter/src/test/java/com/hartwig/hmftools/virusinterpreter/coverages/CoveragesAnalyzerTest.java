package com.hartwig.hmftools.virusinterpreter.coverages;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

import com.google.common.io.Resources;

import org.junit.Test;

public class CoveragesAnalyzerTest {

    private static final String GENOMIC_DIR = Resources.getResource("genomic").getPath();

    private static final String PURPLE_QC_FILE = GENOMIC_DIR + File.separator + "sample.purple.qc";
    private static final String PURPLE_PURITY_TSV = GENOMIC_DIR + File.separator + "sample.purple.purity.tsv";
    private static final String TUMOR_SAMPLE_WGS_METRICS = GENOMIC_DIR + File.separator + "sample.wgsmetrics";

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCalculateExpectedClonalCoverage() throws IOException {
        double expectedClonalCoverage =
                CoveragesAnalyzer.calculateExpectedClonalCoverage(PURPLE_PURITY_TSV, PURPLE_QC_FILE, TUMOR_SAMPLE_WGS_METRICS);
        assertEquals(34.524514945161286, expectedClonalCoverage, EPSILON);
    }
}