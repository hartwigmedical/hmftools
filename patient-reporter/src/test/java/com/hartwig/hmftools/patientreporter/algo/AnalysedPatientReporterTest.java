package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.SampleMetadata;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
import org.junit.Test;

public class AnalysedPatientReporterTest {

    private static final String TUMOR_SAMPLE_ID = "sample";
    private static final String REF_SAMPLE_ID = "ref_sample";
    private static final String PATIENT_ID = "patient";

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_DRIVER_CATALOG_SOMATIC_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_DRIVER_CATALOG_GERMLINE_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.germline.vcf";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_VIRAL_INSERTION_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String PROTECT_EVIDENCE_TSV = BASE_DIRECTORY + "/protect/sample.protect.tsv";
    private static final String PIPELINE_VERSION = BASE_DIRECTORY + "/pipeline.version";
    private static final String MOLECULAR_TISSUE_ORIGIN_PLOT = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = BASE_DIRECTORY + "/cuppa/sample.cuppa.conclusion.txt";
    private static final String VIRAL_BREAKEND_TSV = BASE_DIRECTORY + "/viralbreakend/sample.virusbreakend.vcf.summary.tsv";

    @Test
    public void canRunOnRunDirectory() throws IOException {
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(PatientReporterTestFactory.loadTestAnalysedReportData());

        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .patientId(PATIENT_ID)
                .refSampleId(REF_SAMPLE_ID)
                .refSampleBarcode(Strings.EMPTY)
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .tumorSampleBarcode(Strings.EMPTY)
                .build();

        assertNotNull(reporter.run(sampleMetadata,
                PURPLE_PURITY_TSV,
                PURPLE_QC_FILE,
                PURPLE_DRIVER_CATALOG_SOMATIC_TSV,
                PURPLE_DRIVER_CATALOG_GERMLINE_TSV,
                PURPLE_SOMATIC_VARIANT_VCF,
                PURPLE_GERMLINE_VARIANT_VCF,
                LINX_FUSIONS_TSV,
                LINX_BREAKEND_TSV,
                LINX_VIRAL_INSERTION_TSV,
                LINX_DRIVERS_TSV,
                CHORD_PREDICTION_TXT,
                CIRCOS_FILE,
                PROTECT_EVIDENCE_TSV,
                null,
                false,
                PIPELINE_VERSION,
                MOLECULAR_TISSUE_ORIGIN_TXT,
                MOLECULAR_TISSUE_ORIGIN_PLOT,
                VIRAL_BREAKEND_TSV));
    }

    @Test
    public void canDetermineForNumber() {
        double purityCorrect = 0.40;
        boolean hasReliablePurityCorrect = true;

        Assert.assertEquals(QsFormNumber.FOR_080.display(),
                AnalysedPatientReporter.determineForNumber(hasReliablePurityCorrect, purityCorrect));

        double purityNotCorrect = 0.10;
        boolean hasReliablePurityNotCorrect = false;

        assertEquals(QsFormNumber.FOR_209.display(),
                AnalysedPatientReporter.determineForNumber(hasReliablePurityNotCorrect, purityNotCorrect));
    }
}
