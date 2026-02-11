package com.hartwig.hmftools.orange;

import java.time.LocalDate;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestOrangeConfigFactory
{
    private static final String MELANOMA_DOID = "8923";

    private static final String REFERENCE_SAMPLE_ID = "ref_sample";
    private static final String TUMOR_SAMPLE_ID = "tumor_sample";

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();
    private static final String COHORT_MAPPING_TSV = Resources.getResource("cohort/mapping/example_cohort_mapping.tsv").getPath();
    private static final String DRIVER_GENE_PANEL_TSV = Resources.getResource("driver/example.DriverGenePanel.tsv").getPath();
    private static final String SIGNATURES_ETIOLOGY_TSV =
            Resources.getResource("test_run_resources/sigs/signatures_etiology.tsv").getPath();
    private static final String ENSEMBL_DATA_DIRECTORY = Resources.getResource("test_run_resources/ensembl").getPath();

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String REF_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.summary.tsv";
    private static final String REF_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.flag_counts.tsv";
    private static final String TUMOR_SAMPLE_WGS_METRICS_FILE =
            RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.bam_metric.summary.tsv";
    private static final String TUMOR_SAMPLE_FLAGSTAT_FILE =
            RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.bam_metric.flag_counts.tsv";
    private static final String GERMLINE_GENE_COVERAGE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.gene_coverage.tsv";
    private static final String REF_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/redux/ref_sample.redux.bqr.png";
    private static final String TUMOR_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/redux/tumor_sample.redux.bqr.png";
    private static final String PURPLE_DATA_DIRECTORY = RUN_DIRECTORY + "/purple";
    private static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    private static final String LINX_SOMATIC_DATA_DIRECTORY = RUN_DIRECTORY + "/linx";
    private static final String LINX_GERMLINE_DATA_DIRECTORY = RUN_DIRECTORY + "/linx_germline";
    private static final String LINX_PLOT_DIRECTORY = RUN_DIRECTORY + "/linx/plot";
    private static final String ISOFOX_SUMMARY_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.summary.csv";
    private static final String ISOFOX_GENE_DATA_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.gene_data.csv";
    private static final String ISOFOX_FUSION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.pass_fusions.csv";
    private static final String ISOFOX_ALT_SPLICE_JUNCTION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.alt_splice_junc.csv";
    private static final String LILAC_RESULT_TSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.tsv";
    private static final String LILAC_QC_TSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.qc.tsv";
    private static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusinterprtr/tumor_sample.virus.annotated.tsv";
    private static final String CHORD_PREDICTION_TSV = RUN_DIRECTORY + "/chord/tumor_sample.chord.prediction.tsv";
    private static final String CUPPA_VIS_DATA_TSV = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.vis_data.tsv";
    private static final String CUPPA_SUMMARY_PLOT = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.vis.png";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/ref_sample.peach.haplotypes.best.tsv";
    private static final String SIGS_ALLOCATION_TSV = RUN_DIRECTORY + "/sigs/tumor_sample.sig.allocation.tsv";

    @NotNull
    public static OrangeConfig createMinimalConfig()
    {
        return ImmutableOrangeConfig.builder()
                .experimentType(ExperimentType.TARGETED)
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .samplingDate(LocalDate.now())
                .refGenomeVersion(RefGenomeVersion.V37)
                .outputDir(Strings.EMPTY)
                .doidJsonFile(DOID_JSON)
                .cohortMappingTsv(COHORT_MAPPING_TSV)
                .driverGenePanelTsv(DRIVER_GENE_PANEL_TSV)
                .signaturesEtiologyTsv(SIGNATURES_ETIOLOGY_TSV)
                .ensemblDataDirectory(ENSEMBL_DATA_DIRECTORY)
                .tumorSampleWGSMetricsFile(TUMOR_SAMPLE_WGS_METRICS_FILE)
                .tumorSampleFlagstatFile(TUMOR_SAMPLE_FLAGSTAT_FILE)
                .tumorSampleBqrPlot(TUMOR_SAMPLE_BQR_PLOT)
                .purpleDataDirectory(PURPLE_DATA_DIRECTORY)
                .purplePlotDirectory(PURPLE_PLOT_DIRECTORY)
                .linxSomaticDataDirectory(LINX_SOMATIC_DATA_DIRECTORY)
                .lilacResultTsv(LILAC_RESULT_TSV)
                .lilacQcTsv(LILAC_QC_TSV)
                .limitJsonOutput(false)
                .addDisclaimer(false)
                .build();
    }

    @NotNull
    public static OrangeConfig createTargetedConfig()
    {
        return ImmutableOrangeConfig.builder()
                .from(createMinimalConfig())
                .experimentType(ExperimentType.TARGETED)
                .addPrimaryTumorDoids(MELANOMA_DOID)
                .linxPlotDirectory(LINX_PLOT_DIRECTORY)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .addDisclaimer(true)
                .build();
    }

    @NotNull
    public static OrangeConfig createWGSConfigTumorOnly()
    {
        return ImmutableOrangeConfig.builder()
                .from(createTargetedConfig())
                .experimentType(ExperimentType.WHOLE_GENOME)
                .wgsRefConfig(ImmutableOrangeWGSRefConfig.builder()
                        .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                        .chordPredictionTxt(CHORD_PREDICTION_TSV)
                        .cuppaVisDataTsv(CUPPA_VIS_DATA_TSV)
                        .cuppaSummaryPlot(CUPPA_SUMMARY_PLOT)
                        .sigsAllocationTsv(SIGS_ALLOCATION_TSV)
                        .build())
                .build();
    }

    @NotNull
    public static OrangeConfig createWGSConfigTumorNormal()
    {
        OrangeConfig wgsConfigTumorOnly = createWGSConfigTumorOnly();
        return ImmutableOrangeConfig.builder()
                .from(wgsConfigTumorOnly)
                .rnaConfig(null)
                .wgsRefConfig(ImmutableOrangeWGSRefConfig.builder()
                        .from(wgsConfigTumorOnly.wgsRefConfig())
                        .referenceSampleId(REFERENCE_SAMPLE_ID)
                        .refSampleWGSMetricsFile(REF_SAMPLE_WGS_METRICS_FILE)
                        .refSampleFlagstatFile(REF_SAMPLE_FLAGSTAT_FILE)
                        .germlineGeneCoverageTsv(GERMLINE_GENE_COVERAGE)
                        .refSampleBqrPlot(REF_SAMPLE_BQR_PLOT)
                        .linxGermlineDataDirectory(LINX_GERMLINE_DATA_DIRECTORY)
                        .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                        .build())
                .build();
    }

    @NotNull
    public static OrangeConfig createWGTSConfigTumorNormal()
    {
        // We use tumor_sample as rnaSampleId since we have no real ISOFOX test data for our test_run
        return ImmutableOrangeConfig.builder()
                .from(createWGSConfigTumorNormal())
                .rnaConfig(ImmutableOrangeRnaConfig.builder()
                        .rnaSampleId("tumor_sample")
                        .isofoxSummaryCsv(ISOFOX_SUMMARY_CSV)
                        .isofoxGeneDataCsv(ISOFOX_GENE_DATA_CSV)
                        .isofoxFusionCsv(ISOFOX_FUSION_CSV)
                        .isofoxAltSpliceJunctionCsv(ISOFOX_ALT_SPLICE_JUNCTION_CSV)
                        .build())
                .build();
    }
}
