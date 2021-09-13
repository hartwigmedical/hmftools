package com.hartwig.hmftools.virusinterpreter.coverages;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;

import org.junit.Test;

public class CoveragesAnalyzerTest {

    private static final String GENOMIC_DIR = Resources.getResource("genomic").getPath();

    private static final String PURPLE_QC_FILE = GENOMIC_DIR + File.separator + "sample.purple.qc";
    private static final String PURPLE_PURITY_TSV = GENOMIC_DIR + File.separator + "sample.purple.purity.tsv";
    private static final String TUMOR_SAMPLE_WGS_METRICS = GENOMIC_DIR + File.separator + "sample.wgsmetrics";

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canCalculateExpectedClonalCoverage() throws IOException {
        PurityContext purityContext = PurityContextFile.readWithQC(PURPLE_QC_FILE, PURPLE_PURITY_TSV);

        double expectedClonalCoverage =
                CoveragesAnalyzer.calculateExpectedClonalCoverage(purityContext, TUMOR_SAMPLE_WGS_METRICS);
        assertEquals(34.524514945161286, expectedClonalCoverage, EPSILON);
    }
}