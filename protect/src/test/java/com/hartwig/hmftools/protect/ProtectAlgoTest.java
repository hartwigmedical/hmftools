package com.hartwig.hmftools.protect;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ProtectAlgoTest {

    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline/germline_reporting.tsv").getPath();
    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private static final String SERVE_DIR = Resources.getResource("serve").getPath();

    private static final String PURPLE_PURITY_TSV = Resources.getResource("test_run/purple/sample.purple.purity.tsv").getPath();
    private static final String PURPLE_QC_FILE = Resources.getResource("test_run/purple/sample.purple.qc").getPath();
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.somatic.tsv").getPath();
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.germline.tsv").getPath();
    private static final String PURPLE_SOMATIC_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.somatic.vcf").getPath();
    private static final String PURPLE_GERMLINE_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.germline.vcf").getPath();
    private static final String LINX_FUSION_TSV = Resources.getResource("test_run/linx/sample.linx.fusion.tsv").getPath();
    private static final String LINX_BREAKEND_TSV = Resources.getResource("test_run/linx/sample.linx.breakend.tsv").getPath();
    private static final String LINX_VIRAL_INSERTION_TSV = Resources.getResource("test_run/linx/sample.linx.viral_inserts.tsv").getPath();
    private static final String LINX_DRIVER_CATALOG_TSV = Resources.getResource("test_run/linx/sample.drivers.catalog.tsv").getPath();
    private static final String CHORD_PREDICTION_TXT = Resources.getResource("test_run/chord/sample_chord_prediction.txt").getPath();

    @Test
    public void canRunProtectAlgo() throws IOException {
        ProtectConfig config = ImmutableProtectConfig.builder()
                .tumorSampleId("sample")
                .outputDir(Strings.EMPTY)
                .serveActionabilityDir(SERVE_DIR)
                .doidJsonFile(DOID_JSON)
                .germlineReportingTsv(GERMLINE_REPORTING_TSV)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .linxFusionTsv(LINX_FUSION_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxViralInsertionTsv(LINX_VIRAL_INSERTION_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .build();

        ActionableEvents events = ActionableEventsLoader.readFromDir(SERVE_DIR, RefGenomeVersion.V37);
        ProtectAlgo algo = ProtectAlgo.buildAlgoFromServeActionability(events, Sets.newHashSet());

        assertNotNull(algo.run(config));
    }
}