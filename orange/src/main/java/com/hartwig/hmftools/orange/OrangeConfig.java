package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_FILE_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.util.PathUtil.mandatoryPath;
import static com.hartwig.hmftools.orange.util.PathUtil.optionalPath;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.orange.util.PathResolver;

public class OrangeConfig
{
    public final ExperimentType RunType;

    public final String TumorId;
    public final String ReferenceId;
    public final String RnaSampleId;

    public final RefGenomeVersion RefGenVersion;

    public final Set<String> PrimaryTumorDoids;
    public final LocalDate SamplingDate;

    public final String OutputDir;

    public final String DoidJsonFile;
    public final String SignaturesEtiologyTsv;
    public final String DriverGenePanelTsv;

    public final String PipelineVersionFile;

    public final String PurpleDataDirectory;
    public final String PurplePlotDirectory;

    public final String LinxSomaticDataDirectory;
    public final String LinxGermlineDataDirectory;
    public final String LinxPlotDirectory;

    public final String TumorBamMetricsDir;
    public final String ReferenceBamMetricsDir;

    public final String TumorReduxDir;
    public final String ReferenceReduxDir;

    public final String LilacDir;
    public final String ChordDir;
    public final String CuppaDir;
    public final String PeachDir;
    public final String SigsDir;
    public final String VirusDir;

    public final String IsofoxDir;

    public final boolean LimitJsonOutput;
    public final boolean AddDisclaimer;

    private static final String DOID_SEPARATOR = ";";

    // General params needed for every analysis
    private static final String EXPERIMENT_TYPE = "experiment_type";
    private static final String PRIMARY_TUMOR_DOIDS = "primary_tumor_doids";
    private static final String SAMPLING_DATE = "sampling_date";

    // Input files used by the algorithm
    private static final String DOID_JSON = "doid_json";
    private static final String SIGNATURES_ETIOLOGY_TSV = "signatures_etiology_tsv";

    // Files containing the actual genomic results for this sample.
    private static final String PIPELINE_VERSION_FILE = "pipeline_version_file";

    private static String RNA_SAMPLE_ID = "rna_sample_id";

    private static final String TUMOR_REDUX_DIR_CFG = "tumor_redux_dir";
    private static final String TUMOR_REDUX_DIR_DESC = "Path to Redux tumor files";

    private static String REFERENCE_REDUX_DIR_CFG = "ref_redux_dir";
    private static String REFERENCE_REDUX_DIR_DESC = "Path to Redux reference files";

    // Some additional optional params and flags
    private static final String LIMIT_JSON_OUTPUT = "limit_json_output";
    private static final String ADD_DISCLAIMER = "add_disclaimer";

    public OrangeConfig(final ConfigBuilder configBuilder)
    {
        TumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);

        PrimaryTumorDoids = Sets.newHashSet();

        if(configBuilder.hasValue(PRIMARY_TUMOR_DOIDS))
        {
            String[] values = configBuilder.getValue(PRIMARY_TUMOR_DOIDS).split(DOID_SEPARATOR, -1);
            Arrays.stream(values).forEach(x -> PrimaryTumorDoids.add(x));
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        DoidJsonFile = configBuilder.getValue(DOID_JSON);
        DriverGenePanelTsv = configBuilder.getValue(DRIVER_GENE_PANEL);
        SignaturesEtiologyTsv = configBuilder.getValue(SIGNATURES_ETIOLOGY_TSV);
        PipelineVersionFile = configBuilder.getValue(PIPELINE_VERSION_FILE);

        OutputDir = parseMandatoryOutputDir(configBuilder);

        PathResolver pathResolver = new PathResolver(configBuilder,
                configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR),
                configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        PipelineToolDirectories defaultToolDirectories = resolveDefaultPipelineDirectories(configBuilder);

        PurpleDataDirectory = pathResolver.resolveMandatoryToolDirectory(PURPLE_DIR_CFG, defaultToolDirectories.purpleDir());
        PurplePlotDirectory = pathResolver.resolveMandatoryToolPlotsDirectory(PURPLE_PLOT_DIR_CFG, defaultToolDirectories.purpleDir());
        LinxSomaticDataDirectory = pathResolver.resolveMandatoryToolDirectory(LINX_DIR_CFG, defaultToolDirectories.linxSomaticDir());
        LinxPlotDirectory = optionalPath(pathResolver.resolveOptionalToolPlotsDirectory(LINX_PLOT_DIR_CFG, defaultToolDirectories.linxSomaticDir()));

        if(ReferenceId != null)
        {
            LinxGermlineDataDirectory = pathResolver.resolveMandatoryToolDirectory(
                    LINX_GERMLINE_DIR_CFG, defaultToolDirectories.linxGermlineDir());

            ReferenceBamMetricsDir = pathResolver.resolveMandatoryToolDirectory(REF_METRICS_DIR_CFG, defaultToolDirectories.germlineMetricsDir());

            ReferenceReduxDir = configBuilder.getValue(REFERENCE_REDUX_DIR_CFG);
        }
        else
        {
            LinxGermlineDataDirectory = null;
            ReferenceBamMetricsDir = null;
            ReferenceReduxDir = null;
        }

        TumorBamMetricsDir = pathResolver.resolveMandatoryToolDirectory(TUMOR_METRICS_DIR_CFG, defaultToolDirectories.tumorMetricsDir());
        TumorReduxDir = configBuilder.getValue(TUMOR_REDUX_DIR_CFG);

        LilacDir = pathResolver.resolveOptionalToolDirectory(LILAC_DIR_CFG, defaultToolDirectories.lilacDir());
        ChordDir = pathResolver.resolveMandatoryToolDirectory(CHORD_DIR_CFG, defaultToolDirectories.chordDir());
        CuppaDir = pathResolver.resolveOptionalToolDirectory(CUPPA_DIR_CFG, defaultToolDirectories.cuppaDir());
        PeachDir = pathResolver.resolveOptionalToolDirectory(PEACH_DIR_CFG, defaultToolDirectories.peachDir());
        SigsDir = pathResolver.resolveMandatoryToolDirectory(SIGS_DIR_CFG, defaultToolDirectories.sigsDir());
        VirusDir = pathResolver.resolveOptionalToolDirectory(VIRUS_DIR_CFG, defaultToolDirectories.virusInterpreterDir());

        if(!configBuilder.hasValue(RNA_SAMPLE_ID) || !configBuilder.hasValue(ISOFOX_DIR_CFG))
        {
            RnaSampleId = null;
            IsofoxDir = null;

            LOGGER.info("RNA config not present, will continue without RNA configuration");
        }
        else
        {
            RnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
            IsofoxDir = pathResolver.resolveMandatoryToolDirectory(ISOFOX_DIR_CFG, defaultToolDirectories.isofoxDir());
        }

        LOGGER.debug("RNA sample configured as {}", RnaSampleId);

        AddDisclaimer = configBuilder.hasFlag(ADD_DISCLAIMER);
        if(AddDisclaimer)
        {
            LOGGER.debug("disclaimer will be included in footer");
        }

        LimitJsonOutput = configBuilder.hasFlag(LIMIT_JSON_OUTPUT);
        if(LimitJsonOutput)
        {
            LOGGER.info("JSON limitation has been enabled");
        }

        if(configBuilder.hasValue(SAMPLING_DATE))
        {
            SamplingDate = interpretSamplingDateParam(configBuilder.getValue(SAMPLING_DATE));
        }
        else
        {
            LOGGER.debug("No sampling date has been configured. Setting sampling data to current date");
            SamplingDate = LocalDate.now();
        }

        RunType = determineExperimentType(configBuilder.getValue(EXPERIMENT_TYPE));
        LOGGER.info("experiment type has been resolved to '{}'", RunType);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(EXPERIMENT_TYPE, true, "The type of the experiment, one of WGS or PANEL");
        configBuilder.addConfigItem(TUMOR, true, TUMOR_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addConfigItem(PRIMARY_TUMOR_DOIDS,
                true,
                "A semicolon-separated list of DOIDs representing the primary tumor of patient");
        configBuilder.addConfigItem(SAMPLING_DATE, false, "Optional, if provided represents the sampling date in YYMMDD format");

        addRefGenomeVersion(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.addPath(DOID_JSON, true, "Path to JSON file containing the full DOID tree");
        configBuilder.addPath(SIGNATURES_ETIOLOGY_TSV, true, "Path to signatures etiology TSV");
        addGenePanelOption(configBuilder, true);

        configBuilder.addPath(PIPELINE_VERSION_FILE, false, "Path towards the pipeline version file.");

        // tool output
        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, false, TUMOR_METRICS_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);

        // per tool directory config options are supported, but simpler is to specific the root sample directory containing all tool
        // subdirectories or a single directory containing all pipeline output
        configBuilder.addPath(PIPELINE_SAMPLE_ROOT_DIR, false, PIPELINE_SAMPLE_ROOT_DESC);
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(TUMOR_REDUX_DIR_CFG, false, TUMOR_REDUX_DIR_DESC);
        configBuilder.addPath(REFERENCE_REDUX_DIR_CFG, false, REFERENCE_REDUX_DIR_DESC);

        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_DESC);
        configBuilder.addPath(CHORD_DIR_CFG, false, CHORD_DIR_DESC);
        configBuilder.addPath(CUPPA_DIR_CFG, false, CUPPA_DIR_DESC);
        configBuilder.addPath(PEACH_DIR_CFG, false, PEACH_DIR_DESC);
        configBuilder.addPath(SIGS_DIR_CFG, false, SIGS_DIR_DESC);

        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);
        configBuilder.addPath(SAGE_GERMLINE_DIR_CFG, false, SAGE_GERMLINE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(PURPLE_PLOT_DIR_CFG, false, PURPLE_PLOT_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_PLOT_DIR_CFG, false, LINX_PLOT_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, false, LILAC_DIR_DESC);
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);

        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);

        configBuilder.addFlag(LIMIT_JSON_OUTPUT, "If set, limits every list in the json output to 1 entry.");
        configBuilder.addFlag(ADD_DISCLAIMER, "If set, prints a disclaimer on each page.");
        addLoggingOptions(configBuilder);
    }

    public OrangeRefGenomeVersion orangeRefGenomeVersion()
    {
        return RefGenVersion.is37() ? OrangeRefGenomeVersion.V37 : OrangeRefGenomeVersion.V38;
    }

    private static LocalDate interpretSamplingDateParam(final String samplingDateString)
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

    private static ExperimentType determineExperimentType(final String experimentTypeString)
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

    private static String parseMandatoryOutputDir(final ConfigBuilder configBuilder)
    {
        String dir = parseOutputDir(configBuilder);
        if(dir == null)
        {
            throw new IllegalArgumentException("Could not parse output directory from configuration");
        }
        return mandatoryPath(dir);
    }

    private static PipelineToolDirectories resolveDefaultPipelineDirectories(final ConfigBuilder configBuilder)
    {
        String tumorSampleId = configBuilder.getValue(TUMOR);

        if(configBuilder.hasValue(REFERENCE))
        {
            String referenceSampleId = configBuilder.getValue(REFERENCE);
            return PipelineToolDirectories.resolveToolDirectories(
                    configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG, tumorSampleId, referenceSampleId);
        }
        else
        {
            return PipelineToolDirectories.resolveToolDirectories(
                    configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG, tumorSampleId);
        }
    }

    @VisibleForTesting
    public OrangeConfig(
            final ExperimentType runType, final String tumorId, final String referenceId, final String rnaSampleId,
            final RefGenomeVersion refGenVersion, final Set<String> primaryTumorDoids, final LocalDate samplingDate, final String outputDir,
            final String doidJsonFile, final String signaturesEtiologyTsv, final String driverGenePanelTsv,
            final String pipelineVersionFile, final String purpleDataDirectory, final String purplePlotDirectory,
            final String linxSomaticDataDirectory, final String linxGermlineDataDirectory, final String linxPlotDirectory,
            final String tumorBamMetricsDir, final String referenceBamMetricsDir, final String tumorReduxDir,
            final String referenceReduxDir, final String lilacDir, final String chordDir, final String cuppaDir, final String peachDir,
            final String sigsDir, final String virusDir, final String isofoxDir, final boolean limitJsonOutput, final boolean addDisclaimer)
    {
        RunType = runType;
        TumorId = tumorId;
        ReferenceId = referenceId;
        RnaSampleId = rnaSampleId;
        RefGenVersion = refGenVersion;
        PrimaryTumorDoids = primaryTumorDoids;
        SamplingDate = samplingDate;
        OutputDir = outputDir;
        DoidJsonFile = doidJsonFile;
        SignaturesEtiologyTsv = signaturesEtiologyTsv;
        DriverGenePanelTsv = driverGenePanelTsv;
        PipelineVersionFile = pipelineVersionFile;
        PurpleDataDirectory = purpleDataDirectory;
        PurplePlotDirectory = purplePlotDirectory;
        LinxSomaticDataDirectory = linxSomaticDataDirectory;
        LinxGermlineDataDirectory = linxGermlineDataDirectory;
        LinxPlotDirectory = linxPlotDirectory;
        TumorBamMetricsDir = tumorBamMetricsDir;
        ReferenceBamMetricsDir = referenceBamMetricsDir;
        TumorReduxDir = tumorReduxDir;
        ReferenceReduxDir = referenceReduxDir;
        LilacDir = lilacDir;
        ChordDir = chordDir;
        CuppaDir = cuppaDir;
        PeachDir = peachDir;
        SigsDir = sigsDir;
        VirusDir = virusDir;
        IsofoxDir = isofoxDir;
        LimitJsonOutput = limitJsonOutput;
        AddDisclaimer = addDisclaimer;
    }
}
