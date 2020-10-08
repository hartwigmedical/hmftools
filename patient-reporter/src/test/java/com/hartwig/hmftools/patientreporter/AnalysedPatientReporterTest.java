package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testAnalysedReportData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class AnalysedPatientReporterTest {

    private static final String TUMOR_SAMPLE_ID = "sample";
    private static final String REF_SAMPLE_ID = "ref_sample";

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_GENE_CNV_TSV = BASE_DIRECTORY + "/purple/sample.purple.cnv.gene.tsv";
    private static final String PURPLE_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.tsv";
    private static final String SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String BACHELOR_TSV = BASE_DIRECTORY + "/bachelor/sample.reportable_germline_variants.tsv";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_DISRUPTIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_VIRAL_INSERTION_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.drivers.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";

    @Test
    public void canRunOnRunDirectory() throws IOException {
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(testAnalysedReportData());

        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(REF_SAMPLE_ID)
                .refSampleBarcode(Strings.EMPTY)
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .tumorSampleBarcode(Strings.EMPTY)
                .build();

        assertNotNull(reporter.run(sampleMetadata,
                PURPLE_PURITY_TSV,
                PURPLE_QC_FILE,
                PURPLE_GENE_CNV_TSV,
                PURPLE_DRIVER_CATALOG_TSV,
                SOMATIC_VARIANT_VCF,
                BACHELOR_TSV,
                LINX_FUSIONS_TSV,
                LINX_DISRUPTIONS_TSV,
                LINX_VIRAL_INSERTION_TSV,
                LINX_DRIVERS_TSV,
                CHORD_PREDICTION_TXT,
                CIRCOS_FILE,
                null,
                false));
    }

    @Test
    public void canDetermineForNumber() {
        double purityCorrect = 0.40;
        boolean hasReliablePurityCorrect = true;

        assertEquals(QsFormNumber.FOR_080.display(), AnalysedPatientReporter.determineForNumber(hasReliablePurityCorrect, purityCorrect));

        double purityNotCorrect = 0.10;
        boolean hasReliablePurityNotCorrect = false;

        assertEquals(QsFormNumber.FOR_209.display(),
                AnalysedPatientReporter.determineForNumber(hasReliablePurityNotCorrect, purityNotCorrect));
    }
}
