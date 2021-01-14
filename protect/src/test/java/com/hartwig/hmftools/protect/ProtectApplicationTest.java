package com.hartwig.hmftools.protect;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ProtectApplicationTest {

    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline/germline_reporting.tsv").getPath();
    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private static final String SERVE_DIR = Resources.getResource("serve").getPath();

    private static final String PURPLE_PURITY_TSV = Resources.getResource("test_run/purple/sample.purple.purity.tsv").getPath();
    private static final String PURPLE_QC_FILE = Resources.getResource("test_run/purple/sample.purple.qc").getPath();
    private static final String PURPLE_DRIVER_CATALOG_TSV = Resources.getResource("test_run/purple/sample.driver.catalog.tsv").getPath();
    private static final String PURPLE_SOMATIC_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.somatic.vcf").getPath();
    private static final String BACHELOR_TSV = Resources.getResource("test_run/bachelor/sample.reportable_germline_variants.tsv").getPath();
    private static final String LINX_FUSION_TSV  = Resources.getResource("test_run/linx/sample.linx.fusion.tsv").getPath();
    private static final String LINX_BREAKEND_TSV  = Resources.getResource("test_run/linx/sample.linx.breakend.tsv").getPath();
    private static final String LINX_VIRAL_INSERTION_TSV  = Resources.getResource("test_run/linx/sample.linx.viral_inserts.tsv").getPath();
    private static final String LINX_DRIVERS_TSV  = Resources.getResource("test_run/linx/sample.drivers.catalog.tsv").getPath();
    private static final String CHORD_PREDICTION_TXT  = Resources.getResource("test_run/chord/sample_chord_prediction.txt").getPath();

    @Test
    public void canLoadAndRun() throws IOException {
        ProtectConfig config = ImmutableProtectConfig.builder()
                .tumorSampleId("sample")
                .outputDir(Strings.EMPTY)
                .serveActionabilityDir(SERVE_DIR)
                .doidJsonFile(DOID_JSON)
                .germlineReportingTsv(GERMLINE_REPORTING_TSV)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleDriverCatalogTsv(PURPLE_DRIVER_CATALOG_TSV )
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF )
                .bachelorTsv(BACHELOR_TSV )
                .linxFusionTsv(LINX_FUSION_TSV )
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxViralInsertionTsv(LINX_VIRAL_INSERTION_TSV)
                .linxDriversTsv(LINX_DRIVERS_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .build();

        assertNotNull(ProtectApplication.protectEvidence(config));
    }

}