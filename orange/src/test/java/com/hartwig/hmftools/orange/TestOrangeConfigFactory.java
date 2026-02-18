package com.hartwig.hmftools.orange;

import java.time.LocalDate;
import java.util.Collections;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestOrangeConfigFactory
{
    public static final String REFERENCE_SAMPLE_ID = "ref_sample";
    public static final String TUMOR_SAMPLE_ID = "tumor_sample";
    public static final String RNA_SAMPLE_ID = "rna_sample";

    public static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();
    public static final String MELANOMA_DOID = "8923";
    public static final String COHORT_MAPPING_TSV = Resources.getResource("cohort/mapping/example_cohort_mapping.tsv").getPath();
    public static final String DRIVER_GENE_PANEL_TSV = Resources.getResource("driver/example.DriverGenePanel.tsv").getPath();
    public static final String SIGNATURES_ETIOLOGY_TSV =
            Resources.getResource("test_run_resources/sigs/signatures_etiology.tsv").getPath();

    public static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    public static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    public static final String REF_SAMPLE_WGS_METRICS_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.summary.tsv";
    public static final String REF_SAMPLE_FLAGSTAT_FILE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.flag_counts.tsv";

    public static final String TUMOR_SAMPLE_BAM_METRICS_DIR = RUN_DIRECTORY + "/tumor_sample/bam_metrics/";
    public static final String REFERENCE_SAMPLE_BAM_METRICS_DIR = RUN_DIRECTORY + "/ref_sample/bam_metrics/";
    public static final String TUMOR_SAMPLE_REDUX_DIR = RUN_DIRECTORY + "/redux/";
    public static final String REFERENCE_SAMPLE_REDUX_DIR = RUN_DIRECTORY + "/redux/";

    public static final String TUMOR_SAMPLE_WGS_METRICS_FILE =
            RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.bam_metric.summary.tsv";
    public static final String TUMOR_SAMPLE_FLAGSTAT_FILE =
            RUN_DIRECTORY + "/tumor_sample/bam_metrics/tumor_sample.bam_metric.flag_counts.tsv";
    public static final String GERMLINE_GENE_COVERAGE = RUN_DIRECTORY + "/ref_sample/bam_metrics/ref_sample.bam_metric.gene_coverage.tsv";
    public static final String REF_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/redux/ref_sample.redux.bqr.png";
    public static final String TUMOR_SAMPLE_BQR_PLOT = RUN_DIRECTORY + "/redux/tumor_sample.redux.bqr.png";
    public static final String PURPLE_DATA_DIRECTORY = RUN_DIRECTORY + "/purple";
    public static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    public static final String LINX_SOMATIC_DATA_DIRECTORY = RUN_DIRECTORY + "/linx";
    public static final String LINX_GERMLINE_DATA_DIRECTORY = RUN_DIRECTORY + "/linx_germline";
    public static final String LINX_PLOT_DIRECTORY = RUN_DIRECTORY + "/linx/plot";
    public static final String ISOFOX_DIR = RUN_DIRECTORY + "/isofox/";
    public static final String ISOFOX_SUMMARY_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.isf.summary.tsv";
    public static final String ISOFOX_GENE_DATA_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.isf.gene_data.tsv";
    public static final String ISOFOX_FUSION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.isf.pass_fusions.tsv";
    public static final String ISOFOX_ALT_SPLICE_JUNCTION_CSV = RUN_DIRECTORY + "/isofox/tumor_sample.isf.alt_splice_junc.tsv";
    public static final String LILAC_DIR = RUN_DIRECTORY + "/lilac/";
    public static final String LILAC_RESULT_TSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.tsv";
    public static final String LILAC_QC_TSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.qc.tsv";
    public static final String VIRUS_DIR = RUN_DIRECTORY + "/virusinterprtr/";
    public static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusinterprtr/tumor_sample.virus.annotated.tsv";
    public static final String CHORD_DIR = RUN_DIRECTORY + "/chord/";
    public static final String CHORD_PREDICTION_TSV = RUN_DIRECTORY + "/chord/tumor_sample.chord.prediction.tsv";
    public static final String CUPPA_DIR = RUN_DIRECTORY + "/cuppa/";
    public static final String CUPPA_VIS_DATA_TSV = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.vis_data.tsv";
    public static final String CUPPA_SUMMARY_PLOT = RUN_DIRECTORY + "/cuppa/tumor_sample.cuppa.vis.png";
    public static final String PEACH_DIR = RUN_DIRECTORY + "/peach/";
    public static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/ref_sample.peach.haplotypes.best.tsv";
    public static final String SIGS_DIR = RUN_DIRECTORY + "/sigs/";
    public static final String SIGS_ALLOCATION_TSV = RUN_DIRECTORY + "/sigs/tumor_sample.sig.allocation.tsv";

    public static OrangeConfig createMinimalConfig()
    {
        return new OrangeConfig(
                ExperimentType.TARGETED, TUMOR_SAMPLE_ID, null, null,
                RefGenomeVersion.V37, Collections.emptySet(), LocalDate.now(),
                "", DOID_JSON, COHORT_MAPPING_TSV, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, false);
    }

    public static OrangeConfig createTargetedConfig()
    {
        return new OrangeConfig(
                ExperimentType.TARGETED, TUMOR_SAMPLE_ID, null, null,
                RefGenomeVersion.V37, Set.of(MELANOMA_DOID), LocalDate.now(),
                "", DOID_JSON, COHORT_MAPPING_TSV, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);

        /*
        return ImmutableOrangeConfig.builder()
                .from(createMinimalConfig())
                .experimentType(ExperimentType.TARGETED)
                .addPrimaryTumorDoids(MELANOMA_DOID)
                .linxPlotDirectory(LINX_PLOT_DIRECTORY)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .addDisclaimer(true)
                .build();
        */
    }

    public static OrangeConfig createWGSConfigTumorOnly()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, null, null,
                RefGenomeVersion.V37, Collections.emptySet(), LocalDate.now(),
                "", DOID_JSON, COHORT_MAPPING_TSV, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
        /*
        return ImmutableOrangeConfig.builder()
                .from(createTargetedConfig())
                .experimentType(ExperimentType.WHOLE_GENOME)
                .wgsRefConfig(ImmutableOrangeRefConfig.builder()
                        .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                        .chordPredictionTxt(CHORD_PREDICTION_TSV)
                        .cuppaVisDataTsv(CUPPA_VIS_DATA_TSV)
                        .cuppaSummaryPlot(CUPPA_SUMMARY_PLOT)
                        .sigsAllocationTsv(SIGS_ALLOCATION_TSV)
                        .build())
                .build();
        */
    }

    public static OrangeConfig createWGSConfigTumorNormal()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, REFERENCE_SAMPLE_ID, null,
                RefGenomeVersion.V37, Set.of(MELANOMA_DOID), LocalDate.now(),
                "", DOID_JSON, COHORT_MAPPING_TSV, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);

        /*
        OrangeConfig wgsConfigTumorOnly = createWGSConfigTumorOnly();
        return ImmutableOrangeConfig.builder()
                .from(wgsConfigTumorOnly)
                .rnaConfig(null)
                .wgsRefConfig(ImmutableOrangeRefConfig.builder()
                        .from(wgsConfigTumorOnly.ReferenceConfig)
                        .referenceSampleId(REFERENCE_SAMPLE_ID)
                        .refSampleWGSMetricsFile(REF_SAMPLE_WGS_METRICS_FILE)
                        .refSampleFlagstatFile(REF_SAMPLE_FLAGSTAT_FILE)
                        .germlineGeneCoverageTsv(GERMLINE_GENE_COVERAGE)
                        .refSampleBqrPlot(REF_SAMPLE_BQR_PLOT)
                        .linxGermlineDataDirectory(LINX_GERMLINE_DATA_DIRECTORY)
                        .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                        .build())
                .build();
        */
    }

    @NotNull
    public static OrangeConfig createWGTSConfigTumorNormal()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, REFERENCE_SAMPLE_ID, "tumor_sample",
                RefGenomeVersion.V37, Set.of(MELANOMA_DOID), LocalDate.now(),
                "", DOID_JSON, COHORT_MAPPING_TSV, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);

        /*
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
        */
    }
}
