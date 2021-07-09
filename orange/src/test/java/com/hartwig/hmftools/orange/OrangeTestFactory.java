package com.hartwig.hmftools.orange;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class OrangeTestFactory {

    private static final String MELANOMA_DOID = "8923";

    private static final String REFERENCE_SAMPLE_ID = "ref_sample";
    private static final String TUMOR_SAMPLE_ID = "tumor_sample";

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String PURPLE_PURITY_TSV = RUN_DIRECTORY + "/purple/tumor_sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = RUN_DIRECTORY + "/purple/tumor_sample.purple.qc";
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/tumor_sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/tumor_sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = RUN_DIRECTORY + "/purple/tumor_sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = RUN_DIRECTORY + "/purple/tumor_sample.purple.germline.vcf";
    private static final String PURPLE_GENE_COPY_NUMBER_TSV = RUN_DIRECTORY + "/purple/tumor_sample.purple.cnv.gene.tsv";
    private static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    private static final String LINX_FUSIONS_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.breakend.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/tumor_sample_chord_prediction.txt";
    private static final String CUPPA_CONCLUSION_TXT = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.conclusion.txt";
    private static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusbreakend/tumor_sample.virus.annotated.tsv";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/tumor_sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = RUN_DIRECTORY + "/protect/tumor_sample.protect.tsv";

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private OrangeTestFactory() {
    }

    @NotNull
    public static OrangeConfig createTestOrangeConfig() {
        return ImmutableOrangeConfig.builder()
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .referenceSampleId(REFERENCE_SAMPLE_ID)
                .reportGermline(false)
                .addPrimaryTumorDoids(MELANOMA_DOID)
                .outputDir(Strings.EMPTY)
                .doidJsonFile(DOID_JSON)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleGeneCopyNumberTsv(PURPLE_GENE_COPY_NUMBER_TSV)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .purplePlotDirectory(PURPLE_PLOT_DIRECTORY)
                .linxFusionTsv(LINX_FUSIONS_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .cuppaConclusionTxt(CUPPA_CONCLUSION_TXT)
                .cuppaResultCsv("")
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .protectEvidenceTsv(PROTECT_EVIDENCE_TSV)
                .build();
    }
}
