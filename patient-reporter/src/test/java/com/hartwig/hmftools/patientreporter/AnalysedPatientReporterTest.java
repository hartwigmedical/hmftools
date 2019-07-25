package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testAnalysedReportData;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class AnalysedPatientReporterTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String TUMOR_SAMPLE = "sample";
    private static final String REF_SAMPLE = "ref_sample";

    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity";
    private static final String PURPLE_GENE_CNV_TSV = BASE_DIRECTORY + "/purple/sample.purple.gene.cnv";
    private static final String CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String LINX_FUSIONS_TSV = Resources.getResource("test_run/linx/sample.linx.fusions.tsv").getPath();
    private static final String LINX_DISRUPTIONS_TSV = Resources.getResource("test_run/linx/sample.linx.disruptions.tsv").getPath();
    private static final String BACHELOR_CSV = BASE_DIRECTORY + "/bachelor/sample_germline_variants.csv";
    private static final String CHORD_PREDICTION_FILE = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";

    @Test
    public void canRunOnRunDirectory() throws IOException {
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(testAnalysedReportData());

        assertNotNull(reporter.run(TUMOR_SAMPLE,
                REF_SAMPLE,
                PURPLE_PURITY_TSV,
                PURPLE_GENE_CNV_TSV,
                SOMATIC_VARIANT_VCF,
                LINX_FUSIONS_TSV,
                LINX_DISRUPTIONS_TSV,
                BACHELOR_CSV,
                CHORD_PREDICTION_FILE,
                CIRCOS_FILE,
                null, null));
    }
}
