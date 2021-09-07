package com.hartwig.hmftools.orange;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.orange.report.ImmutableReportConfig;
import com.hartwig.hmftools.orange.report.ReportConfig;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class OrangeTestFactory {

    private static final String MELANOMA_DOID = "8923";

    private static final String REFERENCE_SAMPLE_ID = "ref_sample";
    private static final String TUMOR_SAMPLE_ID = "tumor_sample";

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String REF_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.wgsmetrics";
    private static final String REF_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/ref_sample/flagstat/ref_sample.flagstat";
    private static final String TUMOR_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.wgsmetrics";
    private static final String TUMOR_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/tumor_sample/flagstat/tumor_sample.flagstat";
    private static final String SAGE_GERMLINE_GENE_COVERAGE = RUN_DIRECTORY + "/sage_germline/ref_sample.sage.gene.coverage.tsv";
    private static final String SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/sage_somatic/ref_sample.sage.bqr.png";
    private static final String SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/sage_somatic/tumor_sample.sage.bqr.png";
    private static final String PURPLE_PURITY_TSV = RUN_DIRECTORY + "/purple/tumor_sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = RUN_DIRECTORY + "/purple/tumor_sample.purple.qc";
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/tumor_sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/tumor_sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = RUN_DIRECTORY + "/purple/tumor_sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = RUN_DIRECTORY + "/purple/tumor_sample.purple.germline.vcf";
    private static final String PURPLE_GENE_COPY_NUMBER_TSV = RUN_DIRECTORY + "/purple/tumor_sample.purple.cnv.gene.tsv";
    private static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    private static final String LINX_FUSION_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.breakend.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.driver.catalog.tsv";
    private static final String LINX_DRIVER_TSV = RUN_DIRECTORY + "/linx/tumor_sample.linx.drivers.tsv";
    private static final String LINX_PLOT_DIRECTORY = RUN_DIRECTORY + "/linx/plot";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/tumor_sample_chord_prediction.txt";
    private static final String CUPPA_CONCLUSION_TXT = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.conclusion.txt";
    private static final String CUPPA_RESULT_CSV = RUN_DIRECTORY + "/cuppa/tumor_sample.cup.data.csv";
    private static final String CUPPA_SUMMARY_PLOT = RUN_DIRECTORY + "/cuppa/tumor_sample.cup.report.summary.png";
    private static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusbreakend/tumor_sample.virus.annotated.tsv";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/tumor_sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = RUN_DIRECTORY + "/protect/tumor_sample.protect.tsv";

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private OrangeTestFactory() {
    }

    @NotNull
    public static OrangeConfig createTestOrangeConfig() {
        ReportConfig reportConfig = ImmutableReportConfig.builder().reportGermline(true).maxEvidenceLevel(EvidenceLevel.B).build();

        return ImmutableOrangeConfig.builder()
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .referenceSampleId(REFERENCE_SAMPLE_ID)
                .reportConfig(reportConfig)
                .addPrimaryTumorDoids(MELANOMA_DOID)
                .outputDir(Strings.EMPTY)
                .doidJsonFile(DOID_JSON)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .refSampleWGSMetricsFile(REF_SAMPLE_WGS_METRICS_FILE)
                .refSampleFlagstatFile(REF_SAMPLE_FLAGSTAT_FILE)
                .tumorSampleWGSMetricsFile(TUMOR_SAMPLE_WGS_METRICS_FILE)
                .tumorSampleFlagstatFile(TUMOR_SAMPLE_FLAGSTAT_FILE)
                .sageGermlineGeneCoverageTsv(SAGE_GERMLINE_GENE_COVERAGE)
                .sageSomaticRefSampleBQRPlot(SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT)
                .sageSomaticTumorSampleBQRPlot(SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleGeneCopyNumberTsv(PURPLE_GENE_COPY_NUMBER_TSV)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .purplePlotDirectory(PURPLE_PLOT_DIRECTORY)
                .linxFusionTsv(LINX_FUSION_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .linxDriverTsv(LINX_DRIVER_TSV)
                .linxPlotDirectory(LINX_PLOT_DIRECTORY)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .cuppaConclusionTxt(CUPPA_CONCLUSION_TXT)
                .cuppaResultCsv(CUPPA_RESULT_CSV)
                .cuppaSummaryPlot(CUPPA_SUMMARY_PLOT)
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .protectEvidenceTsv(PROTECT_EVIDENCE_TSV)
                .build();
    }
}
