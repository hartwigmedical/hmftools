package com.hartwig.hmftools.orange;

import java.time.LocalDate;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestOrangeConfigFactory {

    private static final String MELANOMA_DOID = "8923";

    private static final String REFERENCE_SAMPLE_ID = "ref_sample";
    private static final String TUMOR_SAMPLE_ID = "tumor_sample";

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();
    private static final String COHORT_MAPPING_TSV = Resources.getResource("cohort/mapping/example_cohort_mapping.tsv").getPath();
    private static final String COHORT_PERCENTILES_TSV =
            Resources.getResource("cohort/percentile/example_cohort_percentiles.tsv").getPath();
    private static final String DRIVER_GENE_PANEL_TSV = Resources.getResource("driver/example.DriverGenePanel.tsv").getPath();
    private static final String KNOWN_FUSION_FILE = Resources.getResource("known_fusion_data/known_fusion_file.csv").getPath();
    private static final String ISOFOX_GENE_DISTRIBUTION_CSV = Resources.getResource("isofox/empty.gene_distribution.csv").getPath();
    private static final String ISOFOX_ALT_SJ_COHORT_CSV = Resources.getResource("isofox/empty.alt_sj.cohort.csv").getPath();

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String REF_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.wgsmetrics";
    private static final String REF_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/ref_sample/flagstat/ref_sample.flagstat";
    private static final String TUMOR_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.wgsmetrics";
    private static final String TUMOR_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/tumor_sample/flagstat/tumor_sample.flagstat";
    private static final String SAGE_GERMLINE_GENE_COVERAGE = RUN_DIRECTORY + "/sage_germline/ref_sample.sage.gene.coverage.tsv";
    private static final String SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/sage_somatic/ref_sample.sage.bqr.png";
    private static final String SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/sage_somatic/tumor_sample.sage.bqr.png";
    private static final String PURPLE_DATA_DIRECTORY = RUN_DIRECTORY + "/purple";
    private static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    private static final String LINX_SOMATIC_DATA_DIRECTORY = RUN_DIRECTORY + "/linx";
    private static final String LINX_GERMLINE_DATA_DIRECTORY = RUN_DIRECTORY + "/linx_germline";
    private static final String LINX_PLOT_DIRECTORY = RUN_DIRECTORY + "/linx/plot";
    private static final String ISOFOX_SUMMARY_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.summary.csv";
    private static final String ISOFOX_GENE_DATA_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.gene_data.csv";
    private static final String ISOFOX_FUSION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.pass_fusions.csv";
    private static final String ISOFOX_ALT_SPLICE_JUNCTION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.alt_splice_junc.csv";
    private static final String LILAC_RESULT_CSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.csv";
    private static final String LILAC_QC_CSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.qc.csv";
    private static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusbreakend/tumor_sample.virus.annotated.tsv";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/tumor_sample_chord_prediction.txt";
    private static final String CUPPA_RESULT_CSV = RUN_DIRECTORY + "/cuppa/tumor_sample.cup.data.csv";
    private static final String CUPPA_SUMMARY_PLOT = RUN_DIRECTORY + "/cuppa/tumor_sample.cup.report.summary.png";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/tumor_sample.peach.genotype.tsv";

    private TestOrangeConfigFactory() {
    }

    @NotNull
    public static OrangeConfig createDNAConfigTumorOnly() {
        return ImmutableOrangeConfig.builder()
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .addPrimaryTumorDoids(MELANOMA_DOID)
                .experimentDate(LocalDate.now())
                .refGenomeVersion(RefGenomeVersion.V37)
                .outputDir(Strings.EMPTY)
                .doidJsonFile(DOID_JSON)
                .cohortMappingTsv(COHORT_MAPPING_TSV)
                .cohortPercentilesTsv(COHORT_PERCENTILES_TSV)
                .driverGenePanelTsv(DRIVER_GENE_PANEL_TSV)
                .knownFusionFile(KNOWN_FUSION_FILE)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .tumorSampleWGSMetricsFile(TUMOR_SAMPLE_WGS_METRICS_FILE)
                .tumorSampleFlagstatFile(TUMOR_SAMPLE_FLAGSTAT_FILE)
                .sageSomaticTumorSampleBQRPlot(SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT)
                .purpleDataDirectory(PURPLE_DATA_DIRECTORY)
                .purplePlotDirectory(PURPLE_PLOT_DIRECTORY)
                .linxSomaticDataDirectory(LINX_SOMATIC_DATA_DIRECTORY)
                .linxPlotDirectory(LINX_PLOT_DIRECTORY)
                .lilacResultCsv(LILAC_RESULT_CSV)
                .lilacQcCsv(LILAC_QC_CSV)
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .cuppaResultCsv(CUPPA_RESULT_CSV)
                .cuppaSummaryPlot(CUPPA_SUMMARY_PLOT)
                .convertGermlineToSomatic(false)
                .limitJsonOutput(false)
                .build();
    }

    @NotNull
    public static OrangeConfig createDNAConfigTumorNormal() {
        return ImmutableOrangeConfig.builder()
                .from(createDNAConfigTumorOnly())
                .referenceSampleId(REFERENCE_SAMPLE_ID)
                .rnaConfig(null)
                .refSampleWGSMetricsFile(REF_SAMPLE_WGS_METRICS_FILE)
                .refSampleFlagstatFile(REF_SAMPLE_FLAGSTAT_FILE)
                .sageGermlineGeneCoverageTsv(SAGE_GERMLINE_GENE_COVERAGE)
                .sageSomaticRefSampleBQRPlot(SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT)
                .linxGermlineDataDirectory(LINX_GERMLINE_DATA_DIRECTORY)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .build();
    }

    @NotNull
    public static OrangeConfig createDNARNAConfig() {
        // We use tumor_sample as rnaSampleId since we have no real ISOFOX test data for our test_run
        return ImmutableOrangeConfig.builder()
                .from(createDNAConfigTumorNormal())
                .rnaConfig(ImmutableOrangeRNAConfig.builder()
                        .rnaSampleId("tumor_sample")
                        .isofoxGeneDistributionCsv(ISOFOX_GENE_DISTRIBUTION_CSV)
                        .isofoxAltSjCohortCsv(ISOFOX_ALT_SJ_COHORT_CSV)
                        .isofoxSummaryCsv(ISOFOX_SUMMARY_CSV)
                        .isofoxGeneDataCsv(ISOFOX_GENE_DATA_CSV)
                        .isofoxFusionCsv(ISOFOX_FUSION_CSV)
                        .isofoxAltSpliceJunctionCsv(ISOFOX_ALT_SPLICE_JUNCTION_CSV)
                        .build())
                .build();
    }
}
