package com.hartwig.hmftools.orange;

import static java.lang.String.format;

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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.QSEE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.QSEE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_REF_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_REF_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_PLOT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_TUMOR_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_TUMOR_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.util.PathUtil.optionalPath;

import java.nio.file.Files;
import java.nio.file.Paths;
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

    public final String DisplaySampleId;
    public final String PanelName;

    public final RefGenomeVersion RefGenVersion;

    public final Set<String> PrimaryTumorDoids;
    public final String PrimaryTumorLocation;
    public final LocalDate SamplingDate;

    public final String OutputDir;
    public final String OutputId;

    public final String DoidJsonFile;

    public final String PipelineVersionFile;

    public final String PurpleDataDirectory;
    public final String PurplePlotDirectory;
    public final String QSeeDirectory;
    public final String LinxSomaticDataDirectory;
    public final String LinxGermlineDataDirectory;
    public final String LinxPlotDirectory;
    public final String SagePlotDirectory;
    public final String LilacDir;
    public final String ChordDir;
    public final String CuppaDir;
    public final String PeachDir;
    public final String SigsDir;
    public final String VirusDir;
    public final String IsofoxDir;

    public final boolean AddDisclaimer;

    private static final String DOID_SEPARATOR = ";";

    // General params needed for every analysis
    private static final String EXPERIMENT_TYPE = "experiment_type";
    private static final String PRIMARY_TUMOR_DOIDS = "primary_tumor_doids";
    private static final String PRIMARY_TUMOR_LOCATION = "primary_tumor_location";
    private static final String SAMPLING_DATE = "sampling_date";
    private static final String DISPLAY_SAMPLE_ID = "display_sample";
    private static final String PANEL_NAME = "panel_name";

    // Input files used by the algorithm
    private static final String DOID_JSON = "doid_json";

    // Files containing the actual genomic results for this sample.
    private static final String PIPELINE_VERSION_FILE = "pipeline_version_file";

    private static String RNA_SAMPLE_ID = "rna_sample_id";

    // Some additional optional params and flags
    private static final String ADD_DISCLAIMER = "add_disclaimer";

    public OrangeConfig(final ConfigBuilder configBuilder)
    {
        TumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);

        PrimaryTumorDoids = Sets.newHashSet();

        if(configBuilder.hasValue(PRIMARY_TUMOR_DOIDS) && configBuilder.hasValue(DOID_JSON))
        {
            String[] values = configBuilder.getValue(PRIMARY_TUMOR_DOIDS).split(DOID_SEPARATOR, -1);
            Arrays.stream(values).forEach(x -> PrimaryTumorDoids.add(x));
            DoidJsonFile = configBuilder.getValue(DOID_JSON);
        }
        else
        {
            DoidJsonFile = null;
        }

        PrimaryTumorLocation = configBuilder.getValue(PRIMARY_TUMOR_LOCATION);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        String pipelineVersionFile = configBuilder.getValue(PIPELINE_VERSION_FILE);

        String pipelineSampleRootDir = checkAddDirSeparator(configBuilder.getValue(PIPELINE_SAMPLE_ROOT_DIR));

        if(pipelineVersionFile == null && pipelineSampleRootDir != null)
        {
            String testPipelineVersionFile = format("%s/orange_pipeline.version.txt", pipelineSampleRootDir);

            if(Files.exists(Paths.get(testPipelineVersionFile)))
                pipelineVersionFile = testPipelineVersionFile;
        }

        PipelineVersionFile = pipelineVersionFile;

        String outputDir;

        if(pipelineSampleRootDir != null && !configBuilder.hasValue(OUTPUT_DIR))
            outputDir = pipelineSampleRootDir + "orange/";
        else
            outputDir = parseOutputDir(configBuilder);

        OutputDir = outputDir;
        checkCreateOutputDir(OutputDir);

        OutputId = configBuilder.getValue(OUTPUT_ID);

        PathResolver pathResolver = new PathResolver(
                configBuilder,
                pipelineSampleRootDir,
                configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        PipelineToolDirectories defaultToolDirectories = resolveDefaultPipelineDirectories(configBuilder);

        PurpleDataDirectory = pathResolver.resolveMandatoryToolDirectory(PURPLE_DIR_CFG, defaultToolDirectories.purpleDir());
        PurplePlotDirectory = pathResolver.resolveMandatoryToolPlotsDirectory(PURPLE_PLOT_DIR_CFG, defaultToolDirectories.purpleDir());
        LinxSomaticDataDirectory = pathResolver.resolveMandatoryToolDirectory(LINX_DIR_CFG, defaultToolDirectories.linxSomaticDir());
        LinxPlotDirectory = optionalPath(pathResolver.resolveOptionalToolPlotsDirectory(LINX_PLOT_DIR_CFG, defaultToolDirectories.linxSomaticDir()));
        SagePlotDirectory = optionalPath(pathResolver.resolveOptionalToolPlotsDirectory(SAGE_PLOT_DIR_CFG, defaultToolDirectories.linxSomaticDir()));

        QSeeDirectory = pathResolver.resolveOptionalToolDirectory(QSEE_DIR_CFG, defaultToolDirectories.qseeDir());

        if(ReferenceId != null)
        {
            LinxGermlineDataDirectory = pathResolver.resolveMandatoryToolDirectory(
                    LINX_GERMLINE_DIR_CFG, defaultToolDirectories.linxGermlineDir());
        }
        else
        {
            LinxGermlineDataDirectory = null;
        }

        LilacDir = pathResolver.resolveOptionalToolDirectory(LILAC_DIR_CFG, defaultToolDirectories.lilacDir());
        ChordDir = pathResolver.resolveOptionalToolDirectory(CHORD_DIR_CFG, defaultToolDirectories.chordDir());
        CuppaDir = pathResolver.resolveOptionalToolDirectory(CUPPA_DIR_CFG, defaultToolDirectories.cuppaDir());
        PeachDir = pathResolver.resolveOptionalToolDirectory(PEACH_DIR_CFG, defaultToolDirectories.peachDir());
        SigsDir = pathResolver.resolveOptionalToolDirectory(SIGS_DIR_CFG, defaultToolDirectories.sigsDir());
        VirusDir = pathResolver.resolveOptionalToolDirectory(VIRUS_DIR_CFG, defaultToolDirectories.virusInterpreterDir());

        // RNA data can be in the appended Sage VCF and/or from Isofox

        RnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);

        String isofoxDir = null;

        if(configBuilder.hasValue(ISOFOX_DIR_CFG))
        {
            isofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG);
        }
        else
        {
            isofoxDir = pathResolver.resolveOptionalToolDirectory(ISOFOX_DIR_CFG, defaultToolDirectories.isofoxDir());

            if(!Files.exists(Paths.get(isofoxDir)))
                isofoxDir = null;
        }

        IsofoxDir = isofoxDir;

        if(IsofoxDir != null || RnaSampleId != null)
        {
            LOGGER.debug("RNA sample({}) Isofox results({})", RnaSampleId, IsofoxDir);
        }

        DisplaySampleId = configBuilder.getValue(DISPLAY_SAMPLE_ID, TumorId);
        PanelName = configBuilder.getValue(PANEL_NAME);

        AddDisclaimer = configBuilder.hasFlag(ADD_DISCLAIMER);
        if(AddDisclaimer)
        {
            LOGGER.debug("disclaimer will be included in footer");
        }

        if(configBuilder.hasValue(SAMPLING_DATE))
        {
            SamplingDate = interpretSamplingDateParam(configBuilder.getValue(SAMPLING_DATE));
        }
        else
        {
            LOGGER.debug("defaulting sampling data to current date");
            SamplingDate = LocalDate.now();
        }

        RunType = determineExperimentType(configBuilder.getValue(EXPERIMENT_TYPE));
        LOGGER.info("experiment type has been resolved to '{}'", RunType);
    }

    public boolean hasReference() { return ReferenceId != null; }
    public boolean hasRNA() { return RnaSampleId != null || IsofoxDir != null; }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(EXPERIMENT_TYPE, true, "The type of the experiment, one of WGS or PANEL");
        configBuilder.addConfigItem(TUMOR, true, TUMOR_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addConfigItem(PRIMARY_TUMOR_LOCATION, false, "Primary tumor location displayed in the report");

        configBuilder.addConfigItem(PRIMARY_TUMOR_DOIDS,
                false, "A semicolon-separated list of DOIDs representing the primary tumor of patient");

        configBuilder.addConfigItem(SAMPLING_DATE, false, "Optional, if provided represents the sampling date in YYMMDD format");

        addRefGenomeVersion(configBuilder);
        addOutputOptions(configBuilder);

        configBuilder.addPath(DOID_JSON, false, "Path to JSON file containing the full DOID tree");

        configBuilder.addPath(PIPELINE_VERSION_FILE, false, "Path towards the pipeline version file.");

        // tool output
        configBuilder.addPath(BAM_METRICS_TUMOR_DIR_CFG, false, BAM_METRICS_TUMOR_DIR_DESC);
        configBuilder.addPath(BAM_METRICS_REF_DIR_CFG, false, BAM_METRICS_REF_DIR_DESC);

        // per tool directory config options are supported, but simpler is to specific the root sample directory containing all tool
        // subdirectories or a single directory containing all pipeline output
        configBuilder.addPath(PIPELINE_SAMPLE_ROOT_DIR, false, PIPELINE_SAMPLE_ROOT_DESC);

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC); // consider deprecating

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
        configBuilder.addPath(SAGE_PLOT_DIR_CFG, false, SAGE_PLOT_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, false, LILAC_DIR_DESC);
        configBuilder.addPath(QSEE_DIR_CFG, false, QSEE_DIR_DESC);
        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);

        configBuilder.addConfigItem(DISPLAY_SAMPLE_ID, false, "Sample ID for display instead of Tumor ID for demo purposes");
        configBuilder.addConfigItem(PANEL_NAME, false, "Panel name if known");
        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);

        configBuilder.addFlag(ADD_DISCLAIMER, "Prints a disclaimer on each page");
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
            final String doidJsonFile, final String pipelineVersionFile, final String purpleDataDirectory, final String purplePlotDirectory,
            final String linxSomaticDataDirectory, final String linxGermlineDataDirectory, final String linxPlotDirectory,
            final String lilacDir, final String chordDir, final String cuppaDir, final String peachDir, final String sigsDir,
            final String virusDir, final String isofoxDir, final boolean addDisclaimer)
    {
        RunType = runType;
        TumorId = tumorId;
        ReferenceId = referenceId;
        RnaSampleId = rnaSampleId;
        RefGenVersion = refGenVersion;
        PrimaryTumorDoids = primaryTumorDoids;
        SamplingDate = samplingDate;
        OutputDir = outputDir;
        OutputId = null;
        DoidJsonFile = doidJsonFile;
        PrimaryTumorLocation = "";
        PipelineVersionFile = pipelineVersionFile;
        PurpleDataDirectory = purpleDataDirectory;
        PurplePlotDirectory = purplePlotDirectory;
        QSeeDirectory = null;
        LinxSomaticDataDirectory = linxSomaticDataDirectory;
        LinxGermlineDataDirectory = linxGermlineDataDirectory;
        LinxPlotDirectory = checkAddDirSeparator(linxPlotDirectory);
        SagePlotDirectory = linxPlotDirectory; // shared with Linx for testing
        LilacDir = lilacDir;
        ChordDir = chordDir;
        CuppaDir = cuppaDir;
        PeachDir = peachDir;
        SigsDir = sigsDir;
        VirusDir = virusDir;
        IsofoxDir = isofoxDir;
        AddDisclaimer = addDisclaimer;
        DisplaySampleId = tumorId;
        PanelName = "";
    }
}
