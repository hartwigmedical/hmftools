package com.hartwig.hmftools.orange;

import java.time.LocalDate;
import java.util.Collections;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;

import org.jetbrains.annotations.NotNull;

public final class TestOrangeConfigFactory
{
    public static final String REFERENCE_SAMPLE_ID = "ref_sample";
    public static final String TUMOR_SAMPLE_ID = "tumor_sample";
    public static final String RNA_SAMPLE_ID = "rna_sample";

    public static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();
    public static final String MELANOMA_DOID = "8923";
    public static final String DRIVER_GENE_PANEL_TSV = Resources.getResource("driver/example.DriverGenePanel.tsv").getPath();
    public static final String SIGNATURES_ETIOLOGY_TSV =
            Resources.getResource("test_run_resources/sigs/signatures_etiology.tsv").getPath();

    public static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    public static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";

    public static final String TUMOR_SAMPLE_BAM_METRICS_DIR = RUN_DIRECTORY + "/tumor_sample/bam_metrics/";
    public static final String REFERENCE_SAMPLE_BAM_METRICS_DIR = RUN_DIRECTORY + "/ref_sample/bam_metrics/";
    public static final String TUMOR_SAMPLE_REDUX_DIR = RUN_DIRECTORY + "/redux/";
    public static final String REFERENCE_SAMPLE_REDUX_DIR = RUN_DIRECTORY + "/redux/";

    public static final String PURPLE_DATA_DIRECTORY = RUN_DIRECTORY + "/purple";
    public static final String PURPLE_PLOT_DIRECTORY = RUN_DIRECTORY + "/purple/plot";
    public static final String LINX_SOMATIC_DATA_DIRECTORY = RUN_DIRECTORY + "/linx";
    public static final String LINX_GERMLINE_DATA_DIRECTORY = RUN_DIRECTORY + "/linx_germline";
    public static final String LINX_PLOT_DIRECTORY = RUN_DIRECTORY + "/linx/plot";
    public static final String ISOFOX_DIR = RUN_DIRECTORY + "/isofox/";
    public static final String LILAC_DIR = RUN_DIRECTORY + "/lilac/";
    public static final String VIRUS_DIR = RUN_DIRECTORY + "/virusinterprtr/";
    public static final String CHORD_DIR = RUN_DIRECTORY + "/chord/";
    public static final String CUPPA_DIR = RUN_DIRECTORY + "/cuppa/";
    public static final String PEACH_DIR = RUN_DIRECTORY + "/peach/";
    public static final String SIGS_DIR = RUN_DIRECTORY + "/sigs/";

    public static OrangeConfig createMinimalConfig()
    {
        return new OrangeConfig(
                ExperimentType.TARGETED, TUMOR_SAMPLE_ID, null, null,
                RefGenomeVersion.V37, Collections.emptySet(), LocalDate.now(),
                "", DOID_JSON, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
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
                "", DOID_JSON, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
    }

    public static OrangeConfig createWGSConfigTumorOnly()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, null, null,
                RefGenomeVersion.V37, Collections.emptySet(), LocalDate.now(),
                "", DOID_JSON, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
    }

    public static OrangeConfig createWGSConfigTumorNormal()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, REFERENCE_SAMPLE_ID, null,
                RefGenomeVersion.V37, Set.of(MELANOMA_DOID), LocalDate.now(),
                "", DOID_JSON, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
    }

    @NotNull
    public static OrangeConfig createWGTSConfigTumorNormal()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, REFERENCE_SAMPLE_ID, "tumor_sample",
                RefGenomeVersion.V37, Set.of(MELANOMA_DOID), LocalDate.now(),
                "", DOID_JSON, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR,SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
    }
}
