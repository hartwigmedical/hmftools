package com.hartwig.hmftools.protect.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.protect.ImmutableProtectConfig;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ProtectAlgoTest {

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private static final String SERVE_DIR = Resources.getResource("serve").getPath();

    private static final String PURPLE_PURITY_TSV = Resources.getResource("test_run/purple/sample.purple.purity.tsv").getPath();
    private static final String PURPLE_QC_FILE = Resources.getResource("test_run/purple/sample.purple.qc").getPath();
    private static final String PURPLE_GENE_COPY_NUMBER_TSV = Resources.getResource("test_run/purple/sample.purple.cnv.gene.tsv").getPath();
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.somatic.tsv").getPath();
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.germline.tsv").getPath();
    private static final String PURPLE_SOMATIC_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.somatic.vcf").getPath();
    private static final String PURPLE_GERMLINE_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.germline.vcf").getPath();

    private static final String LINX_FUSION_TSV = Resources.getResource("test_run/linx/sample.linx.fusion.tsv").getPath();
    private static final String LINX_BREAKEND_TSV = Resources.getResource("test_run/linx/sample.linx.breakend.tsv").getPath();
    private static final String LINX_DRIVER_CATALOG_TSV = Resources.getResource("test_run/linx/sample.linx.driver.catalog.tsv").getPath();
    private static final String ANNOTATED_VIRUS_TSV = Resources.getResource("test_run/virusbreakend/sample.virus.annotated.tsv").getPath();
    private static final String CHORD_PREDICTION_TXT = Resources.getResource("test_run/chord/sample_chord_prediction.txt").getPath();
    private static final String LILAC_RESULT_CSV = Resources.getResource("test_run/lilac/tumor_sample.lilac.csv").getPath();
    private static final String LILAC_QC_CSV = Resources.getResource("test_run/lilac/tumor_sample.lilac.qc.csv").getPath();

    @Test
    public void canRunProtectAlgo() throws IOException {
        ProtectConfig config = ImmutableProtectConfig.builder()
                .tumorSampleId("sample")
                .referenceSampleId("reference")
                .outputDir(Strings.EMPTY)
                .serveActionabilityDir(SERVE_DIR)
                .refGenomeVersion(RefGenomeVersion.V37)
                .doidJsonFile(DOID_JSON)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleGeneCopyNumberTsv(PURPLE_GENE_COPY_NUMBER_TSV)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .linxFusionTsv(LINX_FUSION_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .lilacResultCsv(LILAC_RESULT_CSV)
                .lilacQcCsv(LILAC_QC_CSV)
                .build();

        ActionableEvents events = ActionableEventsLoader.readFromDir(config.serveActionabilityDir(), config.refGenomeVersion());

        ProtectAlgo algo = ProtectAlgo.build(events, Sets.newHashSet("162"));

        assertNotNull(algo.run(config));
    }
}