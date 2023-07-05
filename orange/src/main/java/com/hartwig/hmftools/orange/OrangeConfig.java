package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Locale;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.sage.SageCommon;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.common.virus.VirusBreakendFile;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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

    // Input files used by the algorithm
    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";
    String COHORT_PERCENTILES_TSV = "cohort_percentiles_tsv";

    // Files containing the actual genomic results for this sample.
    String PIPELINE_VERSION_FILE = "pipeline_version_file";
    String REF_SAMPLE_WGS_METRICS_FILE = "ref_sample_wgs_metrics_file";
    String REF_SAMPLE_FLAGSTAT_FILE = "ref_sample_flagstat_file";
    String TUMOR_SAMPLE_WGS_METRICS_FILE = "tumor_sample_wgs_metrics_file";
    String TUMOR_SAMPLE_FLAGSTAT_FILE = "tumor_sample_flagstat_file";
    String ANNOTATED_VIRUS_TSV = "annotated_virus_tsv";
    String CUPPA_CHART_PLOT = "cuppa_chart_plot";
    String PEACH_GENOTYPE_TSV = "peach_genotype_tsv";

    // Some additional optional params and flags
    String CONVERT_GERMLINE_TO_SOMATIC = "convert_germline_to_somatic";
    String LIMIT_JSON_OUTPUT = "limit_json_output";
    String ADD_DISCLAIMER = "add_disclaimer";
    
    // now using standard config options
    /*
    String LOG_DEBUG = "log_debug";
    String DRIVER_GENE_PANEL_TSV = "driver_gene_panel_tsv";
    String KNOWN_FUSION_FILE = "known_fusion_file";
    String ENSEMBL_DATA_DIRECTORY = "ensembl_data_directory";
    String REF_GENOME_VERSION = "ref_genome_version";
    String OUTPUT_DIRECTORY = "output_dir";
    String PURPLE_DATA_DIRECTORY = "purple_data_directory";
    String LINX_SOMATIC_DATA_DIRECTORY = "linx_somatic_data_directory";
    String LINX_GERMLINE_DATA_DIRECTORY = "linx_germline_data_directory";
    String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    String CUPPA_RESULT_CSV = "cuppa_result_csv";
    String CUPPA_SUMMARY_PLOT = "cuppa_summary_plot";
    String CUPPA_FEATURE_PLOT = "cuppa_feature_plot";
    String SIGS_ALLOCATION_TSV = "sigs_allocation_tsv";
    String SAGE_GERMLINE_GENE_COVERAGE_TSV = "sage_germline_gene_coverage_tsv";
    String SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT = "sage_somatic_ref_sample_bqr_plot";
    String SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT = "sage_somatic_tumor_sample_bqr_plot";
    String PURPLE_PLOT_DIRECTORY = "purple_plot_directory";
    String LINX_PLOT_DIRECTORY = "linx_plot_directory";
    */

    @NotNull
    static void registerConfig(final ConfigBuilder configBuilder) {

        configBuilder.addConfigItem(TUMOR_SAMPLE_ID, true, "The sample ID for which ORANGE will run.");
        configBuilder.addConfigItem(REFERENCE_SAMPLE_ID, false, "(Optional) The reference sample of the tumor sample for which ORANGE will run.");
        configBuilder.addConfigItem(PRIMARY_TUMOR_DOIDS, true, "A semicolon-separated list of DOIDs representing the primary tumor of patient.");
        configBuilder.addConfigItem(EXPERIMENT_DATE, false, "Optional, if provided represents the experiment date in YYMMDD format.");

        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);

        configBuilder.addConfigItem(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        configBuilder.addConfigItem(COHORT_MAPPING_TSV, true, "Path to cohort mapping TSV.");
        configBuilder.addConfigItem(COHORT_PERCENTILES_TSV, true, "Path to cohort percentiles TSV.");
        configBuilder.addPath(DRIVER_GENE_PANEL_OPTION, false, DRIVER_GENE_PANEL_OPTION_DESC);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);

        configBuilder.addConfigItem(PIPELINE_VERSION_FILE, true, "Path towards the pipeline version file.");

        // tool output
        configBuilder.addPath(REF_SAMPLE_WGS_METRICS_FILE, true, "Path towards the ref sample WGS metrics file.");
        configBuilder.addPath(REF_SAMPLE_FLAGSTAT_FILE, true, "Path towards the ref sample flagstat file.");
        configBuilder.addPath(TUMOR_SAMPLE_WGS_METRICS_FILE, true, "Path towards the tumor sample WGS metrics file.");
        configBuilder.addPath(TUMOR_SAMPLE_FLAGSTAT_FILE, true, "Path towards the tumor sample flagstat file.");

        configBuilder.addPath(CHORD_DIR_CFG, true, CHORD_DIR_DESC);
        configBuilder.addPath(CUPPA_DIR_CFG, true, CUPPA_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, true, LILAC_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, true, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, true, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(LINX_PLOT_DIR_CFG, true, LINX_PLOT_DIR_DESC);
        configBuilder.addPath(PEACH_DIR_CFG, true, PEACH_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(PURPLE_PLOT_DIR_CFG, true, PURPLE_PLOT_DIR_DESC);
        configBuilder.addPath(SAGE_DIR_CFG, true, SAGE_DIR_DESC);
        configBuilder.addPath(SAGE_GERMLINE_DIR_CFG, true, SAGE_GERMLINE_DIR_DESC);
        configBuilder.addPath(SIGS_DIR_CFG, true, SIGS_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, true, VIRUS_DIR_DESC);

        /*
        configBuilder.addPath(ANNOTATED_VIRUS_TSV, true, "Path towards the annotated virus TSV.");
        configBuilder.addConfigItem(LILAC_RESULT_CSV, true, "Path towards the LILAC result CSV.");
        configBuilder.addConfigItem(LILAC_QC_CSV, true, "Path towards the LILAC QC CSV.");
        configBuilder.addConfigItem(SAGE_GERMLINE_GENE_COVERAGE_TSV, true, "Path towards the SAGE germline gene coverage TSV.");
        configBuilder.addConfigItem(SAGE_SOMATIC_REF_SAMPLE_BQR_PLOT, true, "Path towards the SAGE somatic ref sample BQR plot.");
        configBuilder.addConfigItem(SAGE_SOMATIC_TUMOR_SAMPLE_BQR_PLOT, true, "Path towards the SAGE somatic tumor sample BQR plot.");
        configBuilder.addConfigItem(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT.");
        configBuilder.addPath(CUPPA_CHART_PLOT, true, "Path towards the Cuppa chart plot PNG.");
        configBuilder.addConfigItem(CUPPA_RESULT_CSV, true, "Path towards the Cuppa result CSV.");
        configBuilder.addConfigItem(CUPPA_SUMMARY_PLOT, true, "Path towards the Cuppa report summary plot PNG.");
        configBuilder.addConfigItem(CUPPA_FEATURE_PLOT, true, "Path towards the Cuppa report feature plot PNG, if present.");
        configBuilder.addConfigItem(PEACH_GENOTYPE_TSV, true, "Path towards the peach genotype TSV.");
        configBuilder.addConfigItem(SIGS_ALLOCATION_TSV, true, "Path towards the signatures allocation TSV.");
        */

        configBuilder.addFlag(CONVERT_GERMLINE_TO_SOMATIC, "If set, germline events are converted to somatic events.");
        configBuilder.addFlag(LIMIT_JSON_OUTPUT, "If set, limits every list in the json output to 1 entry.");
        configBuilder.addFlag(ADD_DISCLAIMER, "If set, prints a disclaimer on each page.");
        addLoggingOptions(configBuilder);

        OrangeRNAConfig.registerConfig(configBuilder);
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
    OrangeRefGenomeVersion refGenomeVersion();

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

    @NotNull
    String ensemblDataDirectory();

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

    @Nullable
    String annotatedVirusTsv();

    @Nullable
    String chordPredictionTxt();

    @Nullable
    String cuppaResultCsv();

    @Nullable
    String cuppaSummaryPlot();

    @Nullable
    String cuppaFeaturePlot();

    @Nullable
    String cuppaChartPlot();

    @Nullable
    String peachGenotypeTsv();

    @Nullable
    String sigsAllocationTsv();

    boolean convertGermlineToSomatic();

    boolean limitJsonOutput();

    boolean addDisclaimer();

    @NotNull
    static OrangeConfig createConfig(final ConfigBuilder configBuilder) throws ParseException {

        setLogLevel(configBuilder);

        if(LOGGER.isDebugEnabled())
            LOGGER.debug("Switched root level logging to DEBUG");

        boolean addDisclaimer = configBuilder.hasFlag(ADD_DISCLAIMER);
        if (addDisclaimer) {
            LOGGER.info("Disclaimer will be included in footer.");
        }

        boolean limitJsonOutput = configBuilder.hasFlag(LIMIT_JSON_OUTPUT);
        if (limitJsonOutput) {
            LOGGER.info("JSON limitation has been enabled.");
        }

        boolean convertGermlineToSomatic = configBuilder.hasFlag(CONVERT_GERMLINE_TO_SOMATIC);
        if (convertGermlineToSomatic) {
            LOGGER.info("Germline conversion to somatic has been enabled.");
        }

        String refSampleId = configBuilder.getValue(REFERENCE_SAMPLE_ID);
        if (refSampleId != null) {
            LOGGER.debug("Ref sample has been configured as {}.", refSampleId);
        }

        LocalDate experimentDate;
        if (configBuilder.hasValue(EXPERIMENT_DATE)) {
            experimentDate = interpretExperimentDateParam(configBuilder.getValue(EXPERIMENT_DATE));
        } else {
            experimentDate = LocalDate.now();
        }

        String tumorSampleId = configBuilder.getValue(TUMOR_SAMPLE_ID);

        String sageSomaticDir = configBuilder.getValue(SAGE_DIR_CFG);
        String sageGermlineDir = configBuilder.getValue(SAGE_GERMLINE_DIR_CFG);
        String sageGemlineGeneCoverage = refSampleId != null ? SageCommon.generateGeneCoverageFilename(sageGermlineDir, refSampleId) : null;
        String sageRefSampleBqrPlot = refSampleId != null ? SageCommon.generateBqrPlotFilename(sageSomaticDir, refSampleId) : null;
        String sageTumorSampleBqrPlot = SageCommon.generateBqrPlotFilename(sageSomaticDir, tumorSampleId);

        String lilacDir = configBuilder.getValue(LILAC_DIR_CFG);
        String lilacCoverage = LilacAllele.generateFilename(lilacDir, tumorSampleId);
        String lilacQc = LilacQcData.generateFilename(lilacDir, tumorSampleId);

        String chordDir = configBuilder.getValue(CHORD_DIR_CFG);
        String chordPredictions = ChordDataFile.generateFilename(chordDir, tumorSampleId);

        String cuppaDir = configBuilder.getValue(CUPPA_DIR_CFG);
        String cuppaDataFile = CuppaDataFile.generateFilename(cuppaDir, tumorSampleId);
        String cuppaSummaryPlot = CuppaDataFile.generateReportSummaryPlotFilename(cuppaDir, tumorSampleId);
        String cuppaFeaturesPlot = CuppaDataFile.generateReportFeaturesPlotFilename(cuppaDir, tumorSampleId);
        String cuppaChartPlot = CuppaDataFile.generateChartPlotFilename(cuppaDir, tumorSampleId);

        String sigsDir = configBuilder.getValue(SIGS_DIR_CFG);
        String sigAllocations = SignatureAllocationFile.generateFilename(sigsDir, tumorSampleId);

        String virusDir = configBuilder.getValue(VIRUS_DIR_CFG);
        String virusAnnotations = AnnotatedVirusFile.generateFileName(virusDir, tumorSampleId);

        String peachDir = configBuilder.getValue(PEACH_DIR_CFG);
        String peachGenotype = checkAddDirSeparator(peachDir) + tumorSampleId + ".peach.genotype.tsv";

        OrangeRefGenomeVersion orangeRefGenomeVersion = OrangeRefGenomeVersion.valueOf(RefGenomeVersion.from(configBuilder).name());

        return ImmutableOrangeConfig.builder()
                .tumorSampleId(configBuilder.getValue(TUMOR_SAMPLE_ID))
                .referenceSampleId(refSampleId)
                .rnaConfig(OrangeRNAConfig.createConfig(configBuilder))
                .primaryTumorDoids(toStringSet(configBuilder.getValue(PRIMARY_TUMOR_DOIDS), DOID_SEPARATOR))
                .experimentDate(experimentDate)
                .refGenomeVersion(orangeRefGenomeVersion)
                .outputDir(parseOutputDir(configBuilder))
                .doidJsonFile(configBuilder.getValue(DOID_JSON))
                .cohortMappingTsv(configBuilder.getValue(COHORT_MAPPING_TSV))
                .cohortPercentilesTsv(configBuilder.getValue(COHORT_PERCENTILES_TSV))
                .driverGenePanelTsv(configBuilder.getValue(DRIVER_GENE_PANEL_OPTION))
                .knownFusionFile(configBuilder.getValue(KNOWN_FUSIONS_FILE))
                .ensemblDataDirectory(configBuilder.getValue(ENSEMBL_DATA_DIR))
                .pipelineVersionFile(configBuilder.getValue(PIPELINE_VERSION_FILE))
                .refSampleWGSMetricsFile(configBuilder.getValue(REF_SAMPLE_WGS_METRICS_FILE))
                .refSampleFlagstatFile(configBuilder.getValue(REF_SAMPLE_FLAGSTAT_FILE))
                .tumorSampleWGSMetricsFile(configBuilder.getValue(TUMOR_SAMPLE_WGS_METRICS_FILE))
                .tumorSampleFlagstatFile(configBuilder.getValue(TUMOR_SAMPLE_FLAGSTAT_FILE))
                .sageGermlineGeneCoverageTsv(sageGemlineGeneCoverage)
                .sageSomaticRefSampleBQRPlot(sageRefSampleBqrPlot)
                .sageSomaticTumorSampleBQRPlot(sageTumorSampleBqrPlot)
                .purpleDataDirectory(configBuilder.getValue(PURPLE_DIR_CFG))
                .purplePlotDirectory(configBuilder.getValue(PURPLE_PLOT_DIR_CFG))
                .linxSomaticDataDirectory(configBuilder.getValue(LINX_DIR_CFG))
                .linxGermlineDataDirectory(configBuilder.getValue(LINX_GERMLINE_DIR_CFG))
                .linxPlotDirectory(configBuilder.getValue(LINX_PLOT_DIR_CFG))
                .lilacResultCsv(lilacCoverage)
                .lilacQcCsv(lilacQc)
                .annotatedVirusTsv(virusAnnotations)
                .chordPredictionTxt(chordPredictions)
                .cuppaResultCsv(cuppaDataFile)
                .cuppaSummaryPlot(cuppaSummaryPlot)
                .cuppaFeaturePlot(cuppaFeaturesPlot)
                .cuppaChartPlot(cuppaChartPlot)
                .peachGenotypeTsv(peachGenotype)
                .sigsAllocationTsv(sigAllocations)
                .convertGermlineToSomatic(convertGermlineToSomatic)
                .limitJsonOutput(limitJsonOutput)
                .addDisclaimer(addDisclaimer)
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
