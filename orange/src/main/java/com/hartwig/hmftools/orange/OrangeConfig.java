package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.CHORD_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.CUPPA_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LILAC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LINX_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LINX_SOMATIC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PEACH_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PURPLE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_SOMATIC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SIGS_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.VIRUS_BREAKEND_DIR;
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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
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

import java.io.File;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Locale;
import java.util.Optional;
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
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.orange.util.InputFiles;

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

    // Some additional optional params and flags
    String CONVERT_GERMLINE_TO_SOMATIC = "convert_germline_to_somatic";
    String LIMIT_JSON_OUTPUT = "limit_json_output";
    String ADD_DISCLAIMER = "add_disclaimer";

    @NotNull
    static void registerConfig(final ConfigBuilder configBuilder) {

        configBuilder.addConfigItem(TUMOR_SAMPLE_ID, true, "The sample ID for which ORANGE will run.");
        configBuilder.addConfigItem(REFERENCE_SAMPLE_ID,
                false,
                "(Optional) The reference sample of the tumor sample for which ORANGE will run.");
        configBuilder.addConfigItem(PRIMARY_TUMOR_DOIDS,
                true,
                "A semicolon-separated list of DOIDs representing the primary tumor of patient.");
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

        // per tool directory config options are supported, but simpler is to specific the root sample directory containing all tool
        // subdirectories or a single directory containing all pipeline output
        configBuilder.addPath(PIPELINE_SAMPLE_ROOT_DIR, false, PIPELINE_SAMPLE_ROOT_DESC);
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(CHORD_DIR_CFG, false, CHORD_DIR_DESC);
        configBuilder.addPath(CUPPA_DIR_CFG, false, CUPPA_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, false, LILAC_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(LINX_PLOT_DIR_CFG, false, LINX_PLOT_DIR_DESC);
        configBuilder.addPath(PEACH_DIR_CFG, false, PEACH_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(PURPLE_PLOT_DIR_CFG, false, PURPLE_PLOT_DIR_DESC);
        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);
        configBuilder.addPath(SAGE_GERMLINE_DIR_CFG, false, SAGE_GERMLINE_DIR_DESC);
        configBuilder.addPath(SIGS_DIR_CFG, false, SIGS_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_DESC);

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
    static OrangeConfig createConfig(final ConfigBuilder configBuilder) {

        setLogLevel(configBuilder);

        if (LOGGER.isDebugEnabled()) {
            LOGGER.debug("Switched root level logging to DEBUG");
        }

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

        String pipelineSampleRootDir = checkAddDirSeparator(configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR));
        String sampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        String sageSomaticDir = InputFiles.directoryExists(getToolDirectory(configBuilder,
                pipelineSampleRootDir,
                sampleDataDir,
                SAGE_DIR_CFG,
                SAGE_SOMATIC_DIR));
        String sageGermlineDir = InputFiles.directoryExists(getToolDirectory(configBuilder,
                pipelineSampleRootDir,
                sampleDataDir,
                SAGE_GERMLINE_DIR_CFG,
                SAGE_GERMLINE_DIR));
        String sageGermlineGeneCoverage =
                refSampleId != null ? SageCommon.generateGeneCoverageFilename(sageGermlineDir, refSampleId) : null;
        String sageSomaticRefSampleBqrPlot = refSampleId != null ? SageCommon.generateBqrPlotFilename(sageSomaticDir, refSampleId) : null;
        String sageSomaticTumorSampleBqrPlot = InputFiles.fileExists(SageCommon.generateBqrPlotFilename(sageSomaticDir, tumorSampleId));

        String purpleDir = InputFiles.directoryExists(getToolDirectory(configBuilder,
                pipelineSampleRootDir,
                sampleDataDir,
                PURPLE_DIR_CFG,
                PURPLE_DIR));
        String purplePlotsDir =
                InputFiles.directoryExists(getToolPlotsDirectory(configBuilder, pipelineSampleRootDir, PURPLE_PLOT_DIR_CFG, PURPLE_DIR));
        String linxSomaticDir = InputFiles.directoryExists(getToolDirectory(configBuilder,
                pipelineSampleRootDir,
                sampleDataDir,
                LINX_DIR_CFG,
                LINX_SOMATIC_DIR));
        String linxGermlineDir =
                getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, LINX_GERMLINE_DIR_CFG, LINX_GERMLINE_DIR);
        String linxPlotsDir = InputFiles.directoryExists(getToolPlotsDirectory(configBuilder,
                pipelineSampleRootDir,
                LINX_PLOT_DIR_CFG,
                LINX_SOMATIC_DIR));

        String lilacDir =
                InputFiles.directoryExists(getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, LILAC_DIR_CFG, LILAC_DIR));
        String lilacCoverage = InputFiles.fileExists(LilacAllele.generateFilename(lilacDir, tumorSampleId));
        String lilacQc = InputFiles.fileExists(LilacQcData.generateFilename(lilacDir, tumorSampleId));

        String chordDir =
                InputFiles.directoryExists(getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, CHORD_DIR_CFG, CHORD_DIR));
        String chordPredictions = ChordDataFile.generateFilename(chordDir, tumorSampleId);

        String cuppaDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, CUPPA_DIR_CFG, CUPPA_DIR);
        String cuppaDataFile = CuppaDataFile.generateFilename(cuppaDir, tumorSampleId);
        String cuppaSummaryPlot = CuppaDataFile.generateReportSummaryPlotFilename(cuppaDir, tumorSampleId);
        String cuppaFeaturesPlot = CuppaDataFile.generateReportFeaturesPlotFilename(cuppaDir, tumorSampleId);
        String cuppaChartPlot = CuppaDataFile.generateChartPlotFilename(cuppaDir, tumorSampleId);

        String sigsDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, SIGS_DIR_CFG, SIGS_DIR);
        String sigAllocations = SignatureAllocationFile.generateFilename(sigsDir, tumorSampleId);

        String virusDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, VIRUS_DIR_CFG, VIRUS_BREAKEND_DIR);
        String virusAnnotations = AnnotatedVirusFile.generateFileName(virusDir, tumorSampleId);

        String peachDir = getToolDirectory(configBuilder, pipelineSampleRootDir, sampleDataDir, PEACH_DIR_CFG, PEACH_DIR);
        String peachGenotype = checkAddDirSeparator(peachDir) + tumorSampleId + ".peach.genotype.tsv";

        OrangeRefGenomeVersion orangeRefGenomeVersion = OrangeRefGenomeVersion.valueOf(RefGenomeVersion.from(configBuilder).name());

        String outputDir = Optional.ofNullable(parseOutputDir(configBuilder))
                .orElseThrow(() -> new IllegalArgumentException("Output directory cannot be null"));
        return ImmutableOrangeConfig.builder()
                .tumorSampleId(configBuilder.getValue(TUMOR_SAMPLE_ID))
                .referenceSampleId(refSampleId)
                .rnaConfig(OrangeRNAConfig.createConfig(configBuilder))
                .primaryTumorDoids(toStringSet(configBuilder.getValue(PRIMARY_TUMOR_DOIDS), DOID_SEPARATOR))
                .experimentDate(experimentDate)
                .refGenomeVersion(orangeRefGenomeVersion)
                .outputDir(outputDir)
                .doidJsonFile(InputFiles.fileExists(configBuilder.getValue(DOID_JSON)))
                .cohortMappingTsv(InputFiles.fileExists(configBuilder.getValue(COHORT_MAPPING_TSV)))
                .cohortPercentilesTsv(InputFiles.fileExists(configBuilder.getValue(COHORT_PERCENTILES_TSV)))
                .driverGenePanelTsv(InputFiles.fileExists(configBuilder.getValue(DRIVER_GENE_PANEL_OPTION)))
                .knownFusionFile(InputFiles.fileExists(configBuilder.getValue(KNOWN_FUSIONS_FILE)))
                .ensemblDataDirectory(InputFiles.fileExists(configBuilder.getValue(ENSEMBL_DATA_DIR)))
                .pipelineVersionFile(InputFiles.fileExists(configBuilder.getValue(PIPELINE_VERSION_FILE)))
                .refSampleWGSMetricsFile(configBuilder.getValue(REF_SAMPLE_WGS_METRICS_FILE))
                .refSampleFlagstatFile(configBuilder.getValue(REF_SAMPLE_FLAGSTAT_FILE))
                .tumorSampleWGSMetricsFile(InputFiles.fileExists(configBuilder.getValue(TUMOR_SAMPLE_WGS_METRICS_FILE)))
                .tumorSampleFlagstatFile(InputFiles.fileExists(configBuilder.getValue(TUMOR_SAMPLE_FLAGSTAT_FILE)))
                .sageGermlineGeneCoverageTsv(sageGermlineGeneCoverage)
                .sageSomaticRefSampleBQRPlot(sageSomaticRefSampleBqrPlot)
                .sageSomaticTumorSampleBQRPlot(sageSomaticTumorSampleBqrPlot)
                .purpleDataDirectory(purpleDir)
                .purplePlotDirectory(purplePlotsDir)
                .linxSomaticDataDirectory(linxSomaticDir)
                .linxGermlineDataDirectory(linxGermlineDir)
                .linxPlotDirectory(linxPlotsDir)
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
    static String getToolDirectory(final ConfigBuilder configBuilder, final String pipelineSampleRootDir, final String sampleDataDir,
            final String toolDirConfig, final String pipelineToolDir) {
        if (configBuilder.hasValue(toolDirConfig)) {
            return checkAddDirSeparator(configBuilder.getValue(toolDirConfig));
        }

        if (pipelineSampleRootDir != null) {
            return pipelineSampleRootDir + pipelineToolDir + File.separator;
        }

        return sampleDataDir;
    }

    @NotNull
    static String getToolPlotsDirectory(final ConfigBuilder configBuilder, final String pipelineSampleRootDir, final String toolDirConfig,
            final String pipelineToolDir) {
        if (configBuilder.hasValue(toolDirConfig)) {
            return configBuilder.getValue(toolDirConfig);
        }

        return pipelineSampleRootDir != null ? pipelineSampleRootDir + pipelineToolDir + "/plot/" : null;
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
