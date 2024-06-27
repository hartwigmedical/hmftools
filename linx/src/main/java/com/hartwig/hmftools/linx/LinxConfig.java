package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.loadDriverGenes;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_GERMLINE_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.linx.types.LinxConstants.DEFAULT_PROXIMITY_DISTANCE;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.linx.analysis.AnnotationExtension;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LinxConfig
{
    public final int ProximityDistance;

    public final String SampleDataPath; // if specified, then will be used for output and purple data directory
    public final String OutputDataPath;
    public final String PurpleDataPath;
    public final String SvVcfFile;
    public final RefGenomeVersion RefGenVersion;

    public final String FragileSiteFile;
    public final String LineElementFile;
    public final int ChainingSvLimit; // for analysis and chaining
    public final boolean IsGermline;

    public boolean LogVerbose;
    public final List<AnnotationExtension> AnnotationExtensions;

    public final LinxOutput Output;

    private final List<String> mSampleIds;

    public final List<DriverGene> DriverGenes;
    public final List<String> RestrictedGeneIds; // specific set of genes to process

    public final boolean RunFusions;
    public final boolean RunDrivers;
    public final boolean FailOnMissing;

    public final int Threads;

    public final ConfigBuilder CmdLineConfig;

    // config options
    public static final String VCF_FILE = "sv_vcf";

    // clustering analysis options
    private static final String CLUSTER_BASE_DISTANCE = "proximity_distance";
    private static final String CHAINING_SV_LIMIT = "chaining_sv_limit";
    private static final String ANNOTATION_EXTENSIONS = "annotations";
    private static final String FAIL_ON_MISSING_SAMPLE = "fail_on_missing";

    public static final String GERMLINE = "germline";

    // reference files
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";

    // logging options
    public static final String LOG_VERBOSE = "log_verbose";

    // global Linx logger
    public static final Logger LNX_LOGGER = LogManager.getLogger(LinxConfig.class);

    // TODO: set by config, try to only use from Linx Config or pass where required
    public static RefGenomeVersion REF_GENOME_VERSION = V37;

    public LinxConfig(final ConfigBuilder configBuilder)
    {
        CmdLineConfig = configBuilder;

        mSampleIds = Lists.newArrayList();
        setSamplesFromConfig(configBuilder);

        String svVcfFile = configBuilder.getValue(VCF_FILE, "");

        if(configBuilder.hasValue(SAMPLE_DATA_DIR_CFG))
        {
            SampleDataPath = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));
            PurpleDataPath = SampleDataPath;

            OutputDataPath = isSingleSample() && !configBuilder.hasValue(OUTPUT_DIR) ? SampleDataPath : parseOutputDir(configBuilder);
        }
        else
        {
            SampleDataPath = "";
            PurpleDataPath = configBuilder.getValue(PURPLE_DIR_CFG, "");
            OutputDataPath = parseOutputDir(configBuilder);
        }

        IsGermline = configBuilder.hasFlag(GERMLINE);

        if(svVcfFile.isEmpty())
        {
            if(IsGermline)
            {
                if(mSampleIds.size() == 1)
                    svVcfFile = PurpleCommon.purpleGermlineSvFile(PurpleDataPath, mSampleIds.get(0));
                else
                    svVcfFile = PurpleDataPath + "*" + PURPLE_SV_GERMLINE_VCF_SUFFIX;
            }
            else
            {
                if(mSampleIds.size() == 1)
                    svVcfFile = PurpleCommon.purpleSomaticSvFile(PurpleDataPath, mSampleIds.get(0));
                else
                    svVcfFile = PurpleDataPath + "*" + PURPLE_SV_VCF_SUFFIX;
            }
        }

        SvVcfFile = svVcfFile;

        RunFusions = configBuilder.hasValue(KNOWN_FUSIONS_FILE);

        Output = new LinxOutput(configBuilder, isSingleSample() && !IsGermline);

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        REF_GENOME_VERSION = RefGenVersion; // see TODO above

        ProximityDistance = configBuilder.getInteger(CLUSTER_BASE_DISTANCE);

        FragileSiteFile = configBuilder.getValue(FRAGILE_SITE_FILE);
        LineElementFile = configBuilder.getValue(LINE_ELEMENT_FILE);

        AnnotationExtensions = AnnotationExtension.fromConfig(configBuilder.getValue(ANNOTATION_EXTENSIONS, ""));

        DriverGenes = loadDriverGenes(configBuilder);
        RunDrivers = !DriverGenes.isEmpty();
        FailOnMissing = configBuilder.hasFlag(FAIL_ON_MISSING_SAMPLE);

        LogVerbose = configBuilder.hasFlag(LOG_VERBOSE);
        Threads = parseThreads(configBuilder);

        ChainingSvLimit = configBuilder.getInteger(CHAINING_SV_LIMIT);

        RestrictedGeneIds = Lists.newArrayList();
        if(configBuilder.hasValue(GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));
            LNX_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }
    }

    public final List<String> getSampleIds() { return mSampleIds; }
    public boolean hasMultipleSamples() { return mSampleIds.size() > 1; }
    public boolean isSingleSample() { return mSampleIds.size() == 1; }

    private void setSamplesFromConfig(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
        {
            String configStr = configBuilder.getValue(SAMPLE);

            if(configStr.contains(CSV_DELIM))
                Arrays.stream(configStr.split(CSV_DELIM, -1)).forEach(x -> mSampleIds.add(x));
            else
                mSampleIds.add(configStr);
        }
        else
        {
            mSampleIds.addAll(loadSampleIdsFile(configBuilder));
        }
    }

    public boolean hasValidSampleDataSource(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(VCF_FILE) && !configBuilder.hasValue(SAMPLE_DATA_DIR_CFG) && !configBuilder.hasValue(PURPLE_DIR_CFG))
        {
            LNX_LOGGER.error("missing SV VCF file or sample data directory");
            return false;
        }

        if(!IsGermline)
        {
            if(PurpleDataPath == null || PurpleDataPath.isEmpty())
            {
                LNX_LOGGER.error("missing purple directory");
                return false;
            }
        }
        else
        {
            if(hasMultipleSamples() && !configBuilder.getValue(VCF_FILE).contains("*"))
            {
                LNX_LOGGER.error("VCF filename mask({}) invalid for multiple samples", configBuilder.getValue(VCF_FILE));
                return false;
            }
        }

        return true;
    }

    public LinxConfig(boolean isGermline)
    {
        ProximityDistance = DEFAULT_PROXIMITY_DISTANCE;
        CmdLineConfig = new ConfigBuilder();
        RefGenVersion = V37;
        PurpleDataPath = "";
        OutputDataPath = null;
        SampleDataPath = "";
        SvVcfFile = "";
        IsGermline = isGermline;
        FragileSiteFile = "";
        LineElementFile = "";
        AnnotationExtensions = Lists.newArrayList();
        mSampleIds = Lists.newArrayList();
        LogVerbose = false;
        Output = new LinxOutput();
        ChainingSvLimit = 0;
        DriverGenes = Lists.newArrayList();
        RestrictedGeneIds = Lists.newArrayList();
        RunDrivers = true;
        RunFusions = true;
        FailOnMissing = false;
        Threads = 0;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, "Sample Id, or list separated by ','");
        configBuilder.addConfigItem(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);
        configBuilder.addConfigItem(PURPLE_DIR_CFG, PURPLE_DIR_DESC);
        configBuilder.addConfigItem(SAMPLE_DATA_DIR_CFG, SAMPLE_DATA_DIR_DESC);
        configBuilder.addConfigItem(RefGenomeVersion.REF_GENOME_VERSION, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addConfigItem(VCF_FILE, "Path to the PURPLE structural variant VCF file");
        configBuilder.addPath(DRIVER_GENE_PANEL_OPTION, false, DRIVER_GENE_PANEL_OPTION_DESC);
        configBuilder.addPath(LINE_ELEMENT_FILE, false, "Line elements file");
        configBuilder.addPath(FRAGILE_SITE_FILE, false, "Fragile site file");
        configBuilder.addFlag(GERMLINE, "Process germline SVs");

        configBuilder.addInteger(CLUSTER_BASE_DISTANCE, "Clustering base distance", DEFAULT_PROXIMITY_DISTANCE);
        configBuilder.addInteger(CHAINING_SV_LIMIT, "Max cluster size for chaining", 0);
        configBuilder.addConfigItem(ANNOTATION_EXTENSIONS, "String list of annotations");

        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);

        LinxOutput.addConfig(configBuilder);
        configBuilder.addFlag(FAIL_ON_MISSING_SAMPLE, "Failing all processing in batch mode if any sample is missing");
        configBuilder.addFlag(LOG_VERBOSE, "Log extra detail");
        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
    }
}
