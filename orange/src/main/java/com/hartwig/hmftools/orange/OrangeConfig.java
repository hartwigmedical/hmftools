package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.FLAGSTAT_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LILAC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.LINX_SOMATIC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.METRICS_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PURPLE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.SAGE_SOMATIC_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PIPELINE_SAMPLE_ROOT_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;
import static com.hartwig.hmftools.orange.util.PathUtil.optionalPath;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Locale;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.sage.SageCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.orange.util.PathResolver;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeConfig
{
    String DOID_SEPARATOR = ";";

    // General params needed for every analysis
    String EXPERIMENT_TYPE = "experiment_type";
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String PRIMARY_TUMOR_DOIDS = "primary_tumor_doids";
    String SAMPLING_DATE = "sampling_date";

    // Input files used by the algorithm
    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";
    String COHORT_PERCENTILES_TSV = "cohort_percentiles_tsv";
    String SIGNATURES_ETIOLOGY_TSV = "signatures_etiology_tsv";

    // Files containing the actual genomic results for this sample.
    String PIPELINE_VERSION_FILE = "pipeline_version_file";
    String TUMOR_SAMPLE_WGS_METRICS_FILE = "tumor_sample_wgs_metrics_file";
    String TUMOR_SAMPLE_FLAGSTAT_FILE = "tumor_sample_flagstat_file";

    // Some additional optional params and flags
    String CONVERT_GERMLINE_TO_SOMATIC = "convert_germline_to_somatic";
    String LIMIT_JSON_OUTPUT = "limit_json_output";
    String ADD_DISCLAIMER = "add_disclaimer";

    static void registerConfig(@NotNull ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(EXPERIMENT_TYPE, true, "The type of the experiment, one of WGS or PANEL");
        configBuilder.addConfigItem(TUMOR_SAMPLE_ID, true, "The sample ID for which ORANGE will run.");
        configBuilder.addConfigItem(PRIMARY_TUMOR_DOIDS,
                true,
                "A semicolon-separated list of DOIDs representing the primary tumor of patient.");
        configBuilder.addConfigItem(SAMPLING_DATE, false, "Optional, if provided represents the sampling date in YYMMDD format.");

        addRefGenomeVersion(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.addPath(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        configBuilder.addPath(COHORT_MAPPING_TSV, true, "Path to cohort mapping TSV.");
        configBuilder.addPath(COHORT_PERCENTILES_TSV, true, "Path to cohort percentiles TSV.");
        addGenePanelOption(configBuilder, true);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);

        configBuilder.addPath(PIPELINE_VERSION_FILE, false, "Path towards the pipeline version file.");

        // tool output
        configBuilder.addPath(TUMOR_SAMPLE_WGS_METRICS_FILE, false, "Path towards the tumor sample WGS metrics file.");
        configBuilder.addPath(TUMOR_SAMPLE_FLAGSTAT_FILE, false, "Path towards the tumor sample flagstat file.");

        // per tool directory config options are supported, but simpler is to specific the root sample directory containing all tool
        // subdirectories or a single directory containing all pipeline output
        configBuilder.addPath(PIPELINE_SAMPLE_ROOT_DIR, false, PIPELINE_SAMPLE_ROOT_DESC);
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(PURPLE_PLOT_DIR_CFG, false, PURPLE_PLOT_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_PLOT_DIR_CFG, false, LINX_PLOT_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, false, LILAC_DIR_DESC);

        configBuilder.addFlag(CONVERT_GERMLINE_TO_SOMATIC, "If set, germline events are converted to somatic events.");
        configBuilder.addFlag(LIMIT_JSON_OUTPUT, "If set, limits every list in the json output to 1 entry.");
        configBuilder.addFlag(ADD_DISCLAIMER, "If set, prints a disclaimer on each page.");
        addLoggingOptions(configBuilder);

        OrangeRnaConfig.registerConfig(configBuilder);
        OrangeWGSRefConfig.registerConfig(configBuilder);
    }

    @NotNull
    ExperimentType experimentType();

    @NotNull
    String tumorSampleId();

    @Nullable
    OrangeRnaConfig rnaConfig();

    @Nullable
    OrangeWGSRefConfig wgsRefConfig();

    @NotNull
    Set<String> primaryTumorDoids();

    @NotNull
    LocalDate samplingDate();

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
    String signaturesEtiologyTsv();

    @NotNull
    String knownFusionFile();

    @NotNull
    String ensemblDataDirectory();

    @Nullable
    String pipelineVersionFile();

    @NotNull
    String tumorSampleWGSMetricsFile();

    @NotNull
    String tumorSampleFlagstatFile();

    @NotNull
    String sageSomaticTumorSampleBQRPlot();

    @NotNull
    String purpleDataDirectory();

    @NotNull
    String purplePlotDirectory();

    @NotNull
    String linxSomaticDataDirectory();

    @Nullable
    String linxPlotDirectory();

    @NotNull
    String lilacResultTsv();

    @NotNull
    String lilacQcTsv();

    boolean convertGermlineToSomatic();

    boolean limitJsonOutput();

    boolean addDisclaimer();

    @NotNull
    static OrangeConfig createConfig(@NotNull ConfigBuilder configBuilder)
    {
        setLogLevel(configBuilder);

        if(LOGGER.isDebugEnabled())
        {
            LOGGER.debug("Switched root level logging to DEBUG");
        }

        boolean addDisclaimer = configBuilder.hasFlag(ADD_DISCLAIMER);
        if(addDisclaimer)
        {
            LOGGER.info("Disclaimer will be included in footer.");
        }

        boolean limitJsonOutput = configBuilder.hasFlag(LIMIT_JSON_OUTPUT);
        if(limitJsonOutput)
        {
            LOGGER.info("JSON limitation has been enabled.");
        }

        boolean convertGermlineToSomatic = configBuilder.hasFlag(CONVERT_GERMLINE_TO_SOMATIC);
        if(convertGermlineToSomatic)
        {
            LOGGER.info("Germline conversion to somatic has been enabled.");
        }

        LocalDate samplingDate;
        if(configBuilder.hasValue(SAMPLING_DATE))
        {
            samplingDate = interpretSamplingDateParam(configBuilder.getValue(SAMPLING_DATE));
        }
        else
        {
            LOGGER.debug("No sampling date has been configured. Setting sampling data to current date");
            samplingDate = LocalDate.now();
        }

        ExperimentType experimentType = determineExperimentType(configBuilder.getValue(EXPERIMENT_TYPE));
        LOGGER.info("Experiment type has been resolved to '{}'", experimentType);

        String tumorSampleId = configBuilder.getValue(TUMOR_SAMPLE_ID);

        PathResolver pathResolver = new PathResolver(configBuilder,
                configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR),
                configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        ImmutableOrangeConfig.Builder builder = ImmutableOrangeConfig.builder();

        builder.experimentType(experimentType)
                .tumorSampleId(tumorSampleId)
                .rnaConfig(OrangeRnaConfig.createConfig(configBuilder, pathResolver))
                .primaryTumorDoids(toStringSet(configBuilder.getValue(PRIMARY_TUMOR_DOIDS), DOID_SEPARATOR))
                .samplingDate(samplingDate)
                .refGenomeVersion(OrangeRefGenomeVersion.valueOf(RefGenomeVersion.from(configBuilder).name()))
                .outputDir(parseMandatoryOutputDir(configBuilder))
                .doidJsonFile(configBuilder.getValue(DOID_JSON))
                .cohortMappingTsv(configBuilder.getValue(COHORT_MAPPING_TSV))
                .cohortPercentilesTsv(configBuilder.getValue(COHORT_PERCENTILES_TSV))
                .driverGenePanelTsv(configBuilder.getValue(DRIVER_GENE_PANEL_OPTION))
                .signaturesEtiologyTsv(configBuilder.getValue(SIGNATURES_ETIOLOGY_TSV))
                .knownFusionFile(configBuilder.getValue(KNOWN_FUSIONS_FILE))
                .ensemblDataDirectory(configBuilder.getValue(ENSEMBL_DATA_DIR))
                .pipelineVersionFile(configBuilder.getValue(PIPELINE_VERSION_FILE))
                .tumorSampleWGSMetricsFile(pathResolver.resolveMandatoryMetricsFile(TUMOR_SAMPLE_WGS_METRICS_FILE, METRICS_DIR, tumorSampleId))
                .tumorSampleFlagstatFile(pathResolver.resolveMandatoryMetricsFile(TUMOR_SAMPLE_FLAGSTAT_FILE, FLAGSTAT_DIR, tumorSampleId))
                .purpleDataDirectory(pathResolver.resolveMandatoryToolDirectory(PURPLE_DIR_CFG, PURPLE_DIR))
                .purplePlotDirectory(pathResolver.resolveMandatoryToolPlotsDirectory(PURPLE_PLOT_DIR_CFG, PURPLE_DIR))
                .linxSomaticDataDirectory(pathResolver.resolveMandatoryToolDirectory(LINX_DIR_CFG, LINX_SOMATIC_DIR))
                .linxPlotDirectory(optionalPath(pathResolver.resolveOptionalToolPlotsDirectory(LINX_PLOT_DIR_CFG, LINX_SOMATIC_DIR)))
                .convertGermlineToSomatic(convertGermlineToSomatic)
                .limitJsonOutput(limitJsonOutput)
                .addDisclaimer(addDisclaimer);

        String sageSomaticDir = pathResolver.resolveMandatoryToolDirectory(SAGE_DIR_CFG, SAGE_SOMATIC_DIR);
        builder.sageSomaticTumorSampleBQRPlot(mandatoryPath(SageCommon.generateBqrPlotFilename(sageSomaticDir, tumorSampleId)));

        String lilacDir = pathResolver.resolveMandatoryToolDirectory(LILAC_DIR_CFG, LILAC_DIR);
        builder.lilacResultTsv(mandatoryPath(LilacAllele.generateFilename(lilacDir, tumorSampleId)));
        builder.lilacQcTsv(mandatoryPath(LilacQcData.generateFilename(lilacDir, tumorSampleId)));

        if(experimentType == ExperimentType.WHOLE_GENOME)
        {
            builder.wgsRefConfig(OrangeWGSRefConfig.createConfig(configBuilder, pathResolver));
        }

        return builder.build();
    }

    @NotNull
    static Iterable<String> toStringSet(@NotNull String paramValue, @NotNull String separator)
    {
        return !paramValue.isEmpty() ? Sets.newHashSet(paramValue.split(separator)) : Sets.newHashSet();
    }

    @NotNull
    private static LocalDate interpretSamplingDateParam(@NotNull String samplingDateString)
    {
        String format = "yyMMdd";

        LocalDate samplingDate;
        try
        {
            samplingDate = LocalDate.parse(samplingDateString, DateTimeFormatter.ofPattern(format, Locale.ENGLISH));
            LOGGER.debug("Configured sampling date to {}", samplingDate);
        }
        catch(DateTimeParseException exception)
        {
            samplingDate = LocalDate.now();
            LOGGER.warn("Could not parse configured sampling date '{}'. Expected format is '{}'", samplingDateString, format);
        }
        return samplingDate;
    }

    @NotNull
    private static ExperimentType determineExperimentType(@NotNull String experimentTypeString)
    {
        switch(experimentTypeString)
        {
            case "WGS":
                return ExperimentType.WHOLE_GENOME;
            case "PANEL":
                return ExperimentType.TARGETED;
            default:
                throw new IllegalArgumentException("Invalid experiment type, must be one of WGS, PANEL");
        }
    }

    @NotNull
    private static String parseMandatoryOutputDir(@NotNull ConfigBuilder configBuilder)
    {
        String dir = parseOutputDir(configBuilder);
        if(dir == null)
        {
            throw new IllegalArgumentException("Could not parse output directory from configuration");
        }
        return mandatoryPath(dir);
    }
}
