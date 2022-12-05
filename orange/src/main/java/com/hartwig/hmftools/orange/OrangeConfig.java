package com.hartwig.hmftools.orange;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Locale;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.orange.util.Config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeConfig {

    Logger LOGGER = LogManager.getLogger(OrangeConfig.class);

    String DOID_SEPARATOR = ";";

    // General params needed for every analysis
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String REFERENCE_SAMPLE_ID = "reference_sample_id";
    String PRIMARY_TUMOR_DOIDS = "primary_tumor_doids";
    String EXPERIMENT_DATE = "experiment_date";
    String REF_GENOME_VERSION = "ref_genome_version";
    String OUTPUT_DIRECTORY = "output_dir";

    // Input files used by the algorithm
    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";
    String COHORT_PERCENTILES_TSV = "cohort_percentiles_tsv";
    String DRIVER_GENE_PANEL_TSV = "driver_gene_panel_tsv";
    String KNOWN_FUSION_FILE = "known_fusion_file";

    // Files containing the actual genomic results for this sample.
    String PIPELINE_VERSION_FILE = "pipeline_version_file";
    String REF_SAMPLE_WGS_METRICS_FILE = "ref_sample_wgs_metrics_file";
    String REF_SAMPLE_FLAGSTAT_FILE = "ref_sample_flagstat_file";
    String TUMOR_SAMPLE_WGS_METRICS_FILE = "tumor_sample_wgs_metrics_file";
    String TUMOR_SAMPLE_FLAGSTAT_FILE = "tumor_sample_flagstat_file";
    String SAGE_GERMLINE_GENE_COVERAGE_TSV = "sage_germline_gene_coverage_tsv";
    String SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT = "sage_somatic_ref_sample_bqr_plot";
    String SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT = "sage_somatic_tumor_sample_bqr_plot";
    String PURPLE_DATA_DIRECTORY = "purple_data_directory";
    String PURPLE_PLOT_DIRECTORY = "purple_plot_directory";
    String LINX_SOMATIC_DATA_DIRECTORY = "linx_somatic_data_directory";
    String LINX_GERMLINE_DATA_DIRECTORY = "linx_germline_data_directory";
    String LINX_PLOT_DIRECTORY = "linx_plot_directory";
    String LILAC_RESULT_CSV = "lilac_result_csv";
    String LILAC_QC_CSV = "lilac_qc_csv";
    String ANNOTATED_VIRUS_TSV = "annotated_virus_tsv";
    String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    String CUPPA_RESULT_CSV = "cuppa_result_csv";
    String CUPPA_SUMMARY_PLOT = "cuppa_summary_plot";
    String CUPPA_FEATURE_PLOT = "cuppa_feature_plot";
    String PEACH_GENOTYPE_TSV = "peach_genotype_tsv";

    // Some additional optional params and flags
    String CONVERT_GERMLINE_TO_SOMATIC = "convert_germline_to_somatic";
    String LIMIT_JSON_OUTPUT = "limit_json_output";
    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which ORANGE will run.");
        options.addOption(REFERENCE_SAMPLE_ID, true, "(Optional) The reference sample of the tumor sample for which ORANGE will run.");
        options.addOption(PRIMARY_TUMOR_DOIDS, true, "A semicolon-separated list of DOIDs representing the primary tumor of patient.");
        options.addOption(EXPERIMENT_DATE, true, "Optional, if provided represents the experiment date in YYMMDD format ");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version used in analysis (37 or 38)");
        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the ORANGE output data will be written to.");

        options.addOption(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        options.addOption(COHORT_MAPPING_TSV, true, "Path to cohort mapping TSV.");
        options.addOption(COHORT_PERCENTILES_TSV, true, "Path to cohort percentiles TSV.");
        options.addOption(DRIVER_GENE_PANEL_TSV, true, "Path to the driver gene panel TSV.");
        options.addOption(KNOWN_FUSION_FILE, true, "Path to the known fusion file.");

        options.addOption(PIPELINE_VERSION_FILE, true, "Path towards the pipeline version file.");
        options.addOption(REF_SAMPLE_WGS_METRICS_FILE, true, "Path towards the ref sample WGS metrics file.");
        options.addOption(REF_SAMPLE_FLAGSTAT_FILE, true, "Path towards the ref sample flagstat file.");
        options.addOption(TUMOR_SAMPLE_WGS_METRICS_FILE, true, "Path towards the tumor sample WGS metrics file.");
        options.addOption(TUMOR_SAMPLE_FLAGSTAT_FILE, true, "Path towards the tumor sample flagstat file.");
        options.addOption(SAGE_GERMLINE_GENE_COVERAGE_TSV, true, "Path towards the SAGE germline gene coverage TSV.");
        options.addOption(SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT, true, "Path towards the SAGE somatic ref sample BQR plot.");
        options.addOption(SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT, true, "Path towards the SAGE somatic tumor sample BQR plot.");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Path towards the purple data.");
        options.addOption(PURPLE_PLOT_DIRECTORY, true, "Path towards the directory holding all purple plots.");
        options.addOption(LINX_SOMATIC_DATA_DIRECTORY, true, "Path towards the somatic linx data.");
        options.addOption(LINX_GERMLINE_DATA_DIRECTORY, true, "Path towards the germline linx data.");
        options.addOption(LINX_PLOT_DIRECTORY, true, "Path towards the directory holding all linx plots.");
        options.addOption(LILAC_RESULT_CSV, true, "Path towards the LILAC result CSV.");
        options.addOption(LILAC_QC_CSV, true, "Path towards the LILAC QC CSV.");
        options.addOption(ANNOTATED_VIRUS_TSV, true, "Path towards the annotated virus TSV.");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT.");
        options.addOption(CUPPA_RESULT_CSV, true, "Path towards the Cuppa result CSV.");
        options.addOption(CUPPA_SUMMARY_PLOT, true, "Path towards the Cuppa report summary plot PNG.");
        options.addOption(CUPPA_FEATURE_PLOT, true, "Path towards the Cuppa report feature plot PNG.");
        options.addOption(PEACH_GENOTYPE_TSV, true, "Path towards the peach genotype TSV.");

        options.addOption(CONVERT_GERMLINE_TO_SOMATIC, false, "If set, germline events are converted to somatic events");
        options.addOption(LIMIT_JSON_OUTPUT, false, "If set, limits every list in the json output to 1 entry.");
        options.addOption(LOG_DEBUG, false, "If set, set the log level to debug rather than default.");

        for (Option rnaOption : OrangeRNAConfig.createOptions().getOptions()) {
            options.addOption(rnaOption);
        }

        return options;
    }

    @NotNull
    String tumorSampleId();

    @Nullable
    String referenceSampleId();

    @Nullable
    OrangeRNAConfig rnaConfig();

    @NotNull
    Set<String> primaryTumorDoids();

    @NotNull
    LocalDate experimentDate();

    @NotNull
    RefGenomeVersion refGenomeVersion();

    @NotNull
    String outputDir();

    @NotNull
    String doidJsonFile();

    @NotNull
    String cohortMappingTsv();

    @NotNull
    String cohortPercentilesTsv();

    @NotNull
    String driverGenePanelTsv();

    @NotNull
    String knownFusionFile();

    @Nullable
    String pipelineVersionFile();

    @Nullable
    String refSampleWGSMetricsFile();

    @Nullable
    String refSampleFlagstatFile();

    @NotNull
    String tumorSampleWGSMetricsFile();

    @NotNull
    String tumorSampleFlagstatFile();

    @Nullable
    String sageGermlineGeneCoverageTsv();

    @Nullable
    String sageSomaticRefSampleBQRPlot();

    @NotNull
    String sageSomaticTumorSampleBQRPlot();

    @NotNull
    String purpleDataDirectory();

    @NotNull
    String purplePlotDirectory();

    @NotNull
    String linxSomaticDataDirectory();

    @Nullable
    String linxGermlineDataDirectory();

    @NotNull
    String linxPlotDirectory();

    @NotNull
    String lilacResultCsv();

    @NotNull
    String lilacQcCsv();

    @NotNull
    String annotatedVirusTsv();

    @NotNull
    String chordPredictionTxt();

    @NotNull
    String cuppaResultCsv();

    @NotNull
    String cuppaSummaryPlot();

    @Nullable
    String cuppaFeaturePlot();

    @Nullable
    String peachGenotypeTsv();

    boolean convertGermlineToSomatic();

    boolean limitJsonOutput();

    @NotNull
    static OrangeConfig createConfig(@NotNull CommandLine cmd) throws ParseException, IOException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
            LOGGER.debug("Switched root level logging to DEBUG");
        }

        boolean limitJsonOutput = cmd.hasOption(LIMIT_JSON_OUTPUT);
        if (limitJsonOutput) {
            LOGGER.info("JSON limitation has been enabled.");
        }

        boolean convertGermlineToSomatic = cmd.hasOption(CONVERT_GERMLINE_TO_SOMATIC);
        if (convertGermlineToSomatic) {
            LOGGER.info("Germline events will be converted to somatic events");
        }

        String refSampleId = Config.optionalValue(cmd, REFERENCE_SAMPLE_ID);
        if (refSampleId != null) {
            LOGGER.debug("Ref sample configured as {}", refSampleId);
        }

        LocalDate experimentDate;
        if (cmd.hasOption(EXPERIMENT_DATE)) {
            experimentDate = interpretExperimentDateParam(cmd.getOptionValue(EXPERIMENT_DATE));
        } else {
            experimentDate = LocalDate.now();
        }

        return ImmutableOrangeConfig.builder()
                .tumorSampleId(Config.nonOptionalValue(cmd, TUMOR_SAMPLE_ID))
                .referenceSampleId(refSampleId)
                .rnaConfig(OrangeRNAConfig.createConfig(cmd))
                .primaryTumorDoids(toStringSet(Config.nonOptionalValue(cmd, PRIMARY_TUMOR_DOIDS), DOID_SEPARATOR))
                .experimentDate(experimentDate)
                .refGenomeVersion(RefGenomeVersion.from(Config.nonOptionalValue(cmd, REF_GENOME_VERSION)))
                .outputDir(Config.outputDir(cmd, OUTPUT_DIRECTORY))
                .doidJsonFile(Config.nonOptionalFile(cmd, DOID_JSON))
                .cohortMappingTsv(Config.nonOptionalFile(cmd, COHORT_MAPPING_TSV))
                .cohortPercentilesTsv(Config.nonOptionalFile(cmd, COHORT_PERCENTILES_TSV))
                .driverGenePanelTsv(Config.nonOptionalFile(cmd, DRIVER_GENE_PANEL_TSV))
                .knownFusionFile(Config.nonOptionalFile(cmd, KNOWN_FUSION_FILE))
                .pipelineVersionFile(Config.optionalValue(cmd, PIPELINE_VERSION_FILE))
                .refSampleWGSMetricsFile(Config.optionalFile(cmd, REF_SAMPLE_WGS_METRICS_FILE))
                .refSampleFlagstatFile(Config.optionalFile(cmd, REF_SAMPLE_FLAGSTAT_FILE))
                .tumorSampleWGSMetricsFile(Config.nonOptionalFile(cmd, TUMOR_SAMPLE_WGS_METRICS_FILE))
                .tumorSampleFlagstatFile(Config.nonOptionalFile(cmd, TUMOR_SAMPLE_FLAGSTAT_FILE))
                .sageGermlineGeneCoverageTsv(Config.optionalFile(cmd, SAGE_GERMLINE_GENE_COVERAGE_TSV))
                .sageSomaticRefSampleBQRPlot(Config.optionalFile(cmd, SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT))
                .sageSomaticTumorSampleBQRPlot(Config.nonOptionalFile(cmd, SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT))
                .purpleDataDirectory(Config.nonOptionalDir(cmd, PURPLE_DATA_DIRECTORY))
                .purplePlotDirectory(Config.nonOptionalDir(cmd, PURPLE_PLOT_DIRECTORY))
                .linxSomaticDataDirectory(Config.nonOptionalDir(cmd, LINX_SOMATIC_DATA_DIRECTORY))
                .linxGermlineDataDirectory(Config.optionalDir(cmd, LINX_GERMLINE_DATA_DIRECTORY))
                .linxPlotDirectory(Config.nonOptionalValue(cmd, LINX_PLOT_DIRECTORY))
                .lilacResultCsv(Config.nonOptionalFile(cmd, LILAC_RESULT_CSV))
                .lilacQcCsv(Config.nonOptionalFile(cmd, LILAC_QC_CSV))
                .annotatedVirusTsv(Config.nonOptionalFile(cmd, ANNOTATED_VIRUS_TSV))
                .chordPredictionTxt(Config.nonOptionalFile(cmd, CHORD_PREDICTION_TXT))
                .cuppaResultCsv(Config.nonOptionalFile(cmd, CUPPA_RESULT_CSV))
                .cuppaSummaryPlot(Config.nonOptionalFile(cmd, CUPPA_SUMMARY_PLOT))
                .cuppaFeaturePlot(Config.optionalValue(cmd, CUPPA_FEATURE_PLOT))
                .peachGenotypeTsv(Config.optionalFile(cmd, PEACH_GENOTYPE_TSV))
                .convertGermlineToSomatic(convertGermlineToSomatic)
                .limitJsonOutput(limitJsonOutput)
                .build();
    }

    @NotNull
    static Iterable<String> toStringSet(@NotNull String paramValue, @NotNull String separator) {
        return !paramValue.isEmpty() ? Sets.newHashSet(paramValue.split(separator)) : Sets.newHashSet();
    }

    @NotNull
    private static LocalDate interpretExperimentDateParam(@NotNull String experimentDateString) {
        String format = "yyMMdd";

        LocalDate experimentDate;
        try {
            experimentDate = LocalDate.parse(experimentDateString, DateTimeFormatter.ofPattern(format, Locale.ENGLISH));
            LOGGER.debug("Configured experiment date to {}", experimentDate);
        } catch (DateTimeParseException exception) {
            experimentDate = LocalDate.now();
            LOGGER.warn("Could not parse configured experiment date '{}'. Expected format is '{}'", experimentDateString, format);
        }
        return experimentDate;
    }
}
