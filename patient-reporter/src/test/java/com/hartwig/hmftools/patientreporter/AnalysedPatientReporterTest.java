package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testAnalysedReportData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;
import com.hartwig.hmftools.patientreporter.variants.germline.BachelorFile;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
    private static final String CHORD_PREDICTION_FILE = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String LINX_VIRAL_INSERTIONS_FILE = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final Logger LOGGER = LogManager.getLogger(AnalysedPatientReporterTest.class);

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
                CHORD_PREDICTION_FILE,
                CIRCOS_FILE,
                LINX_VIRAL_INSERTIONS_FILE,
                null,
                false));
    }

    @Test
    public void canMergeViralInsertions() throws IOException {
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(testAnalysedReportData());
        List<ViralInsertion> viralInsertions = reporter.analyzeViralInsertions(LINX_VIRAL_INSERTIONS_FILE);
        assertEquals(2, viralInsertions.size());

        ViralInsertion viralInsertion1 = viralInsertions.get(0);
        assertEquals("Human papillomavirus type 15", viralInsertion1.virus());
        assertEquals("1", viralInsertion1.countVirus());

        ViralInsertion viralInsertion2 = viralInsertions.get(1);
        assertEquals("Human papillomavirus type 16", viralInsertion2.virus());
        assertEquals("2", viralInsertion2.countVirus());

    }
}
