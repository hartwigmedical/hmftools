package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testAnalysedReportData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class AnalysedPatientReporterTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String TUMOR_SAMPLE_ID = "sample";
    private static final String REF_SAMPLE_ID = "ref_sample";

    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity";
    private static final String PURPLE_GENE_CNV_TSV = BASE_DIRECTORY + "/purple/sample.purple.gene.cnv";
    private static final String CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String LINX_FUSIONS_TSV = Resources.getResource("test_run/linx/sample.linx.fusions.tsv").getPath();
    private static final String LINX_DISRUPTIONS_TSV = Resources.getResource("test_run/linx/sample.linx.disruptions.tsv").getPath();
    private static final String BACHELOR_CSV = BASE_DIRECTORY + "/bachelor/sample_germline_variants.csv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String LINX_VIRAL_INSERTIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/linx/sample.drivers.catalog.tsv";

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
                PURPLE_GENE_CNV_TSV,
                SOMATIC_VARIANT_VCF,
                LINX_FUSIONS_TSV,
                LINX_DISRUPTIONS_TSV,
                BACHELOR_CSV,
                CHORD_PREDICTION_TXT,
                CIRCOS_FILE,
                LINX_VIRAL_INSERTIONS_TSV,
                LINX_DRIVER_CATALOG_TSV,
                null,
                false));
    }

}
