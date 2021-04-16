package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class GenomicAnalyzerTest {

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
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.drivers.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String PROTECT_EVIDENCE_TSV = BASE_DIRECTORY + "/protect/sample.protect.tsv";

    @Test
    public void canRunOnTestRun() throws IOException {
        GenomicAnalyzer analyzer = new GenomicAnalyzer();

        assertNotNull(analyzer.run("sample",
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
                PROTECT_EVIDENCE_TSV));
    }
}
