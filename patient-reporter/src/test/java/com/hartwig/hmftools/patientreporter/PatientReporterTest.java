package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSvAnalyzerModel;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;

import org.junit.Test;

public class PatientReporterTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String TUMOR_SAMPLE = "CPCT11111111T";
    private static final String REF_SAMPLE = "CPCT11111111R";

    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/CPCT11111111T.purple.purity";
    private static final String PURPLE_GENE_CNV_TSV = BASE_DIRECTORY + "/purple/CPCT11111111T.purple.gene.cnv";
    private static final String SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/kodu.sage.vcf.gz";
    private static final String BACHELOR_CSV = BASE_DIRECTORY + "/bachelor/CPCT11111111T_germline_variants.csv";
    private static final String CHORD_PREDICTION_FILE = BASE_DIRECTORY + "/chord/CPCT11111111T_chord_prediction.txt";

    private static final String CIRCOS_FILE = Resources.getResource("circos").getPath() + "/circos_example.png";

    @Test
    public void canRunOnRunDirectory() throws IOException {
        final BaseReportData baseReportData = testBaseReportData();
        final SequencedReportData reporterData = testSequencedReportData();
        final SvAnalyzer svAnalyzer = testSvAnalyzerModel();
        final PatientReporter reporter = new PatientReporter(baseReportData, reporterData, svAnalyzer);

        assertNotNull(reporter.run(TUMOR_SAMPLE,
                REF_SAMPLE,
                PURPLE_PURITY_TSV,
                PURPLE_GENE_CNV_TSV,
                SOMATIC_VARIANT_VCF,
                BACHELOR_CSV,
                CHORD_PREDICTION_FILE,
                CIRCOS_FILE,
                null));
    }
}
