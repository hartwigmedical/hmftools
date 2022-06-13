package com.hartwig.hmftools.rose;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class RoseAlgoTest {

    private static final String ACTIONABILITY_DATABASE_TSV = Resources.getResource("actionability/ActionabilityDB.tsv").getPath();

    private static final String PURPLE_PURITY_TSV = Resources.getResource("test_run/purple/sample.purple.purity.tsv").getPath();
    private static final String PURPLE_QC_FILE = Resources.getResource("test_run/purple/sample.purple.qc").getPath();
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.somatic.tsv").getPath();
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV =
            Resources.getResource("test_run/purple/sample.driver.catalog.germline.tsv").getPath();
    private static final String PURPLE_SOMATIC_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.somatic.vcf").getPath();
    private static final String PURPLE_GERMLINE_VARIANT_VCF = Resources.getResource("test_run/purple/sample.purple.germline.vcf").getPath();
    private static final String PURPLE_GENE_COPY_NUMBER_TSV = Resources.getResource("test_run/purple/sample.purple.cnv.gene.tsv").getPath();

    private static final String LINX_FUSION_TSV = Resources.getResource("test_run/linx/sample.linx.fusion.tsv").getPath();
    private static final String LINX_BREAKEND_TSV = Resources.getResource("test_run/linx/sample.linx.breakend.tsv").getPath();
    private static final String LINX_DRIVER_CATALOG_TSV = Resources.getResource("test_run/linx/sample.linx.driver.catalog.tsv").getPath();

    private static final String CHORD_PREDICTION_TXT = Resources.getResource("test_run/chord/sample_chord_prediction.txt").getPath();

    private static final String ANNOTATED_VIRUS_TSV = Resources.getResource("test_run/virusbreakend/sample.virus.annotated.tsv").getPath();

    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = Resources.getResource("test_run/cuppa/sample.cuppa.conclusion.txt").getPath();

    private static final String DRIVER_GENE_37_TSV = Resources.getResource("drivercatalog/driver.gene.panel.tsv").getPath();
    private static final String DRIVER_GENE_38_TSV = Resources.getResource("drivercatalog/driver.gene.panel.tsv").getPath();
    private static final String PRIMARY_TUMOR_TSV = Resources.getResource("primarytumor/primary_tumor.tsv").getPath();

    @Test
    public void canRunProtectAlgo() throws IOException {
        RoseConfig config = ImmutableRoseConfig.builder()
                .outputDir(Strings.EMPTY)
                .actionabilityDatabaseTsv(ACTIONABILITY_DATABASE_TSV)
                .refGenomeVersion(RefGenomeVersion.V37)
                .tumorSampleId("sample")
                .refSampleId("reference")
                .patientId("patientId")
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .purpleSomaticCopyNumberTsv(PURPLE_GENE_COPY_NUMBER_TSV)
                .linxFusionTsv(LINX_FUSION_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .molecularTissueOriginTxt(MOLECULAR_TISSUE_ORIGIN_TXT)
                .driverGeneTsv(DRIVER_GENE_38_TSV)
                .primaryTumorTsv(PRIMARY_TUMOR_TSV)
                .build();

        RoseAlgo algo = RoseAlgo.build(config.actionabilityDatabaseTsv(),
                config.driverGeneTsv(),
                config.refGenomeVersion(),
                config.primaryTumorTsv());

        assertNotNull(algo.run(config));
    }
}