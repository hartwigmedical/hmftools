package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.SvFileLoader.VCF_FILE;
import static com.hartwig.hmftools.linx.germline.GermlineFilter.GERMLINE_MIN_QUAL;
import static com.hartwig.hmftools.linx.types.LinxConstants.DEFAULT_CHAINING_SV_LIMIT;
import static com.hartwig.hmftools.linx.types.LinxConstants.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.hasDatabaseConfig;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.germline.GermlinePonCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LinxConfig
{
    public final int ProximityDistance;

    public final String SampleDataPath; // if specified, then will be used for output and purple data directory
    public final String OutputDataPath;
    public final String PurpleDataPath;
    public final String SvVcfFile;

    public final boolean UploadToDB;
    public final String FragileSiteFile;
    public final String KataegisFile;
    public final String LineElementFile;
    public final int ChainingSvLimit; // for analysis and chaining
    public final boolean IsGermline;
    public final boolean IndelAnnotation;
    public final String IndelFile;

    public boolean LogVerbose;
    public final List<String> RequiredAnnotations;

    public final LinxOutput Output;

    private final List<String> mSampleIds;

    public final List<DriverGene> DriverGenes;
    public final List<String> RestrictedGeneIds; // specific set of genes to process

    public final CommandLine CmdLineArgs;
    public final boolean RunFusions;
    public final boolean RunDrivers;
    public final boolean HomDisAllGenes;

    public final int Threads;

    // config options
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";
    public static final String PURPLE_DATA_DIR = "purple_dir";
    public static final String SAMPLE = "sample";
    public static final String UPLOAD_TO_DB = "upload_to_db"; // true by default when in single-sample mode, false for batch

    public static final String CHECK_DRIVERS = "check_drivers";
    public static final String CHECK_FUSIONS = "check_fusions";
    public static final String HOM_DIS_ALL_GENES = "hom_dis_all_genes";

    // clustering analysis options
    private static final String CLUSTER_BASE_DISTANCE = "proximity_distance";
    private static final String CHAINING_SV_LIMIT = "chaining_sv_limit";
    private static final String REQUIRED_ANNOTATIONS = "annotations";

    public static RefGenomeVersion RG_VERSION = V37;

    private static final String INDEL_ANNOTATIONS = "indel_annotation";
    private static final String GERMLINE = "germline";

    // reference files
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String KATAEGIS_FILE = "kataegis_file";
    private static final String INDEL_FILE = "indel_input_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    public static final String GENE_ID_FILE = "gene_id_file";

    private static final String THREADS = "threads";

    // logging options
    public static final String LOG_VERBOSE = "log_verbose";

    // global Linx logger
    public static final Logger LNX_LOGGER = LogManager.getLogger(LinxConfig.class);

    public LinxConfig(final CommandLine cmd)
    {
        CmdLineArgs = cmd;
        mSampleIds = sampleListFromConfigStr(cmd.getOptionValue(SAMPLE));

        if(cmd.hasOption(UPLOAD_TO_DB) && hasDatabaseConfig(cmd))
        {
            UploadToDB = Boolean.parseBoolean(cmd.getOptionValue(UPLOAD_TO_DB));
        }
        else
        {
            UploadToDB = mSampleIds.size() == 1;
        }

        String svVcfFile = cmd.getOptionValue(VCF_FILE, "");

        if(cmd.hasOption(SAMPLE_DATA_DIR))
        {
            SampleDataPath = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR));
            PurpleDataPath = SampleDataPath;

            OutputDataPath = isSingleSample() ? SampleDataPath : parseOutputDir(cmd);

            if(svVcfFile.isEmpty() && mSampleIds.size() == 1)
            {
                svVcfFile = SampleDataPath + mSampleIds.get(0) + ".purple.sv.vcf.gz";
            }
            else
            {
                svVcfFile = SampleDataPath + "*.purple.sv.vcf.gz";
            }
        }
        else
        {
            SampleDataPath = "";
            PurpleDataPath = cmd.getOptionValue(PURPLE_DATA_DIR, "");
            OutputDataPath = parseOutputDir(cmd);
        }

        SvVcfFile = svVcfFile;

        IsGermline = cmd.hasOption(GERMLINE);

        RunFusions = cmd.hasOption(CHECK_FUSIONS);

        Output = new LinxOutput(cmd, isSingleSample() && !IsGermline);

        if(cmd.hasOption(REF_GENOME_VERSION))
            RG_VERSION = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

        ProximityDistance = cmd.hasOption(CLUSTER_BASE_DISTANCE) ? Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE))
                : DEFAULT_PROXIMITY_DISTANCE;

        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        KataegisFile = cmd.getOptionValue(KATAEGIS_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        IndelAnnotation = cmd.hasOption(INDEL_ANNOTATIONS);
        IndelFile = cmd.getOptionValue(INDEL_FILE, "");

        RequiredAnnotations = Lists.newArrayList();

        if(cmd.hasOption(REQUIRED_ANNOTATIONS))
        {
            Arrays.stream(cmd.getOptionValue(REQUIRED_ANNOTATIONS).split(ITEM_DELIM, -1)).forEach(x -> RequiredAnnotations.add(x));
        }

        DriverGenes = loadDriverGenes(cmd);
        RunDrivers = cmd.hasOption(CHECK_DRIVERS) && !DriverGenes.isEmpty();
        HomDisAllGenes = cmd.hasOption(HOM_DIS_ALL_GENES);

        LogVerbose = cmd.hasOption(LOG_VERBOSE);
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        ChainingSvLimit = cmd.hasOption(CHAINING_SV_LIMIT) ? Integer.parseInt(cmd.getOptionValue(CHAINING_SV_LIMIT)) : DEFAULT_CHAINING_SV_LIMIT;

        RestrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));
            LNX_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }
    }

    public boolean breakendGeneLoading()
    {
        return isSingleSample() && !RunDrivers && RestrictedGeneIds.isEmpty() && !IsGermline;
    }

    private List<DriverGene> loadDriverGenes(final CommandLine cmd)
    {
        if(DriverGenePanelConfig.isConfigured(cmd))
        {
            try
            {
                return DriverGenePanelConfig.driverGenes(cmd);
            }
            catch (IOException e)
            {
                LNX_LOGGER.error("invalid driver gene panel file({})", cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
            }
        }

        return Lists.newArrayList();
    }

    public final List<String> getSampleIds() { return mSampleIds; }
    public void setSampleIds(final List<String> list) { mSampleIds.addAll(list); }

    public boolean hasMultipleSamples() { return mSampleIds.size() > 1; }
    public boolean isSingleSample() { return mSampleIds.size() == 1; }

    public boolean loadSampleDataFromFile() { return (!PurpleDataPath.isEmpty() && SvVcfFile != null) || IsGermline; }

    public static List<String> sampleListFromConfigStr(final String configSampleStr)
    {
        final List<String> sampleIds = Lists.newArrayList();

        if(configSampleStr != null && !configSampleStr.equals("*"))
        {
            if(configSampleStr.contains(","))
            {
                String[] tumorList = configSampleStr.split(",");
                sampleIds.addAll(Arrays.stream(tumorList).collect(Collectors.toList()));
            }
            else if(configSampleStr.contains(".csv"))
            {
                sampleIds.addAll(ConfigUtils.loadSampleIdsFile(configSampleStr));
                LNX_LOGGER.info("Loaded {} specific sample IDs", sampleIds.size());
            }
            else
            {
                // assume refers to a single sample
                sampleIds.add(configSampleStr);
            }
        }

        return sampleIds;
    }

    public boolean hasValidSampleDataSource(final CommandLine cmd)
    {
        if(!cmd.hasOption(VCF_FILE) && !cmd.hasOption(SAMPLE_DATA_DIR))
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
            if(hasMultipleSamples() && !cmd.getOptionValue(VCF_FILE).contains("*"))
            {
                LNX_LOGGER.error("VCF filename mask({}) invalid for multiple samples", cmd.getOptionValue(VCF_FILE));
                return false;
            }
        }

        return true;
    }

    public LinxConfig()
    {
        CmdLineArgs = null;
        ProximityDistance = DEFAULT_PROXIMITY_DISTANCE;
        RG_VERSION = V37;
        PurpleDataPath = "";
        OutputDataPath = null;
        SampleDataPath = "";
        SvVcfFile = "";
        UploadToDB = false;
        IsGermline = false;
        FragileSiteFile = "";
        KataegisFile = "";
        LineElementFile = "";
        IndelFile = "";
        IndelAnnotation = false;
        RequiredAnnotations = Lists.newArrayList();
        mSampleIds = Lists.newArrayList();
        LogVerbose = false;
        Output = new LinxOutput();
        ChainingSvLimit = DEFAULT_CHAINING_SV_LIMIT;
        DriverGenes = Lists.newArrayList();
        RestrictedGeneIds = Lists.newArrayList();
        RunDrivers = false;
        HomDisAllGenes = false;
        RunFusions = false;
        Threads = 0;
    }

    public static boolean validConfig(final CommandLine cmd)
    {
        return configPathValid(cmd, PURPLE_DATA_DIR) && configPathValid(cmd, SAMPLE_DATA_DIR)
            && configPathValid(cmd, FRAGILE_SITE_FILE) && configPathValid(cmd, KATAEGIS_FILE) && configPathValid(cmd, LINE_ELEMENT_FILE)
            && configPathValid(cmd, ENSEMBL_DATA_DIR) && configPathValid(cmd, VCF_FILE) && configPathValid(cmd, DRIVER_GENE_PANEL_OPTION)
            && configPathValid(cmd, INDEL_FILE) && FusionDisruptionAnalyser.validConfig(cmd);
    }

    public static boolean configPathValid(final CommandLine cmd, final String configItem)
    {
        if(!cmd.hasOption(configItem))
            return true;

        final String filePath = cmd.getOptionValue(configItem);

        if(filePath.contains("*")) // too difficult to determine
            return true;

        if(!Files.exists(Paths.get(filePath)))
        {
            LNX_LOGGER.error("invalid config path: {} = {}", configItem, filePath);
            return false;
        }

        return true;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(PURPLE_DATA_DIR, true, "Sample purple data directory");
        addOutputDir(options);
        options.addOption(SAMPLE_DATA_DIR, true, "Optional: directory for per-sample SV data, default is to use output_dir");
        options.addOption(SAMPLE, true, "Sample Id, or list separated by ';' or '*' for all in DB");
        options.addOption(UPLOAD_TO_DB, true, "Upload all LINX data to DB (true/false), single-sample default=true, batch-mode default=false");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default), or 38");
        options.addOption(CHECK_DRIVERS, false, "Check SVs against drivers catalog");
        options.addOption(CHECK_FUSIONS, false, "Run fusion detection");
        options.addOption(HOM_DIS_ALL_GENES, false, "Run fusion detection");
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 5000");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file");
        options.addOption(KATAEGIS_FILE, true, "Kataegis data file");
        options.addOption(GERMLINE, false, "Process germline SVs");
        options.addOption(GERMLINE_MIN_QUAL, true, "Germline SV MinQual limit");
        options.addOption(GENE_ID_FILE, true, "Limit to Ensembl gene ids specified in file");
        options.addOption(CHAINING_SV_LIMIT, true, "Optional: max cluster size for chaining");
        options.addOption(REQUIRED_ANNOTATIONS, true, "Optional: string list of annotations");
        options.addOption(INDEL_ANNOTATIONS, false, "Optional: annotate clusters and TIs with INDELs");
        options.addOption(INDEL_FILE, true, "Optional: cached set of INDELs");
        options.addOption(THREADS, true, "Cohort mode: number of threads");

        ConfigUtils.addLoggingOptions(options);
        options.addOption(LOG_VERBOSE, false, "Log extra detail");

        LinxOutput.addCmdLineArgs(options);
        GermlinePonCache.addCmdLineArgs(options);
    }
}
