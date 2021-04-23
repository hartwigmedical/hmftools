package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.linx.SvFileLoader.VCF_FILE;
import static com.hartwig.hmftools.linx.types.LinxConstants.DEFAULT_CHAINING_SV_LIMIT;
import static com.hartwig.hmftools.linx.types.LinxConstants.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.hasDatabaseConfig;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cli.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;

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
    public final String ReplicationOriginsFile;
    public final String ViralHostsFile;
    public final int ChainingSvLimit; // for analysis and chaining
    public final boolean IsGermline;
    public final boolean IndelAnnotation;
    public final String IndelFile;

    public boolean LogVerbose;
    public String RequiredAnnotations;

    public final LinxOutput Output;

    private final List<String> mSampleIds;

    public final List<DriverGene> DriverGenes;
    public final List<String> RestrictedGeneIds; // specific set of genes to process

    // config options
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";
    public static final String PURPLE_DATA_DIR = "purple_dir";
    public static final String SAMPLE = "sample";
    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String UPLOAD_TO_DB = "upload_to_db"; // true by default when in single-sample mode, false for batch

    public static final String CHECK_DRIVERS = "check_drivers";
    public static final String CHECK_FUSIONS = "check_fusions";

    // clustering analysis options
    private static final String CLUSTER_BASE_DISTANCE = "proximity_distance";
    private static final String CHAINING_SV_LIMIT = "chaining_sv_limit";
    private static final String REQUIRED_ANNOTATIONS = "annotations";

    public static RefGenomeVersion RG_VERSION = V37;

    private static final String INDEL_ANNOTATIONS = "indel_annotation";

    // reference files
    public static final String REF_GENOME_FILE = "ref_genome";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String KATAEGIS_FILE = "kataegis_file";
    private static final String INDEL_FILE = "indel_input_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String VIRAL_HOSTS_FILE = "viral_hosts_file";
    private static final String REPLICATION_ORIGINS_FILE = "replication_origins_file";
    public static final String GENE_ID_FILE = "gene_id_file";
    private static final String GERMLINE = "germline";

    // logging options
    public static final String LOG_DEBUG = "log_debug";
    public static final String LOG_VERBOSE = "log_verbose";

    // global Linx logger
    public static final Logger LNX_LOGGER = LogManager.getLogger(LinxConfig.class);

    public LinxConfig(final CommandLine cmd)
    {
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
            // /DO52664T.purple.sv.vcf.gz
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

        Output = new LinxOutput(cmd, isSingleSample());

        if(cmd.hasOption(REF_GENOME_VERSION))
            RG_VERSION = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

        ProximityDistance = cmd.hasOption(CLUSTER_BASE_DISTANCE) ? Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE))
                : DEFAULT_PROXIMITY_DISTANCE;

        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        KataegisFile = cmd.getOptionValue(KATAEGIS_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        ViralHostsFile = cmd.getOptionValue(VIRAL_HOSTS_FILE, "");
        ReplicationOriginsFile = cmd.getOptionValue(REPLICATION_ORIGINS_FILE, "");
        IndelAnnotation = cmd.hasOption(INDEL_ANNOTATIONS);
        IndelFile = cmd.getOptionValue(INDEL_FILE, "");
        RequiredAnnotations = cmd.getOptionValue(REQUIRED_ANNOTATIONS, "");

        DriverGenes = loadDriverGenes(cmd);

        LogVerbose = cmd.hasOption(LOG_VERBOSE);

        ChainingSvLimit = cmd.hasOption(CHAINING_SV_LIMIT) ? Integer.parseInt(cmd.getOptionValue(CHAINING_SV_LIMIT)) : DEFAULT_CHAINING_SV_LIMIT;

        RestrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            loadGeneIdsFile(inputFile, RestrictedGeneIds);
            LNX_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }
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

    private void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LNX_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            int geneIdIndex = 0;
            boolean parseData = false;
            if (fileContents.get(0).contains("GeneId"))
            {
                final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(fileContents.get(0), ",");
                geneIdIndex = fieldsIndexMap.get("GeneId");
                parseData = fieldsIndexMap.size() > 1;
                fileContents.remove(0);
            }

            for(String data : fileContents)
            {
                if(data.startsWith("#"))
                    continue;

                if(parseData)
                {
                    final String[] items = data.split(",", -1);
                    geneIdList.add(items[geneIdIndex]);
                }
                else
                {
                    geneIdList.add(data);
                }
            }
        }
        catch (IOException e)
        {
            LNX_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
            }
    }


    public final List<String> getSampleIds() { return mSampleIds; }
    public void setSampleIds(final List<String> list) { mSampleIds.addAll(list); }
    public boolean hasMultipleSamples() { return mSampleIds.size() > 1; }
    public boolean isSingleSample() { return mSampleIds.size() == 1; }

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
                sampleIds.addAll(loadSampleListFile(configSampleStr));
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

    public LinxConfig(int proximityDistance)
    {
        ProximityDistance = proximityDistance;
        RG_VERSION = V37;
        PurpleDataPath = "";
        OutputDataPath = "";
        SampleDataPath = "";
        SvVcfFile = "";
        UploadToDB = false;
        IsGermline = false;
        FragileSiteFile = "";
        KataegisFile = "";
        LineElementFile = "";
        ViralHostsFile = "";
        IndelFile = "";
        ReplicationOriginsFile = "";
        IndelAnnotation = false;
        RequiredAnnotations = "";
        mSampleIds = Lists.newArrayList();
        LogVerbose = false;
        Output = new LinxOutput();
        ChainingSvLimit = DEFAULT_CHAINING_SV_LIMIT;
        DriverGenes = Lists.newArrayList();
        RestrictedGeneIds = Lists.newArrayList();
    }

    public static boolean validConfig(final CommandLine cmd)
    {
        return configPathValid(cmd, PURPLE_DATA_DIR) && configPathValid(cmd, SAMPLE_DATA_DIR)
            && configPathValid(cmd, FRAGILE_SITE_FILE) && configPathValid(cmd, KATAEGIS_FILE) && configPathValid(cmd, LINE_ELEMENT_FILE)
            && configPathValid(cmd, GENE_TRANSCRIPTS_DIR) && configPathValid(cmd, VCF_FILE) && configPathValid(cmd, DRIVER_GENE_PANEL_OPTION)
            && configPathValid(cmd, VIRAL_HOSTS_FILE) && configPathValid(cmd, REPLICATION_ORIGINS_FILE) && configPathValid(cmd, INDEL_FILE)
            && FusionDisruptionAnalyser.validConfig(cmd);
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
        options.addOption(OUTPUT_DIR, true, "Linx output directory");
        options.addOption(SAMPLE_DATA_DIR, true, "Optional: directory for per-sample SV data, default is to use output_dir");
        options.addOption(SAMPLE, true, "Sample Id, or list separated by ';' or '*' for all in DB");
        options.addOption(UPLOAD_TO_DB, true, "Upload all LINX data to DB (true/false), single-sample default=true, batch-mode default=false");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default), or 38");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, "Driver gene panel file");
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 5000");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file");
        options.addOption(VIRAL_HOSTS_FILE, true, "Viral hosts file");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file");
        options.addOption(KATAEGIS_FILE, true, "Kataegis data file");
        options.addOption(REPLICATION_ORIGINS_FILE, true, "Origins of replication file");
        options.addOption(GERMLINE, false, "Process germline SVs");
        options.addOption(GENE_ID_FILE, true, "Limit to Ensembl gene ids specified in file");
        options.addOption(CHAINING_SV_LIMIT, true, "Optional: max cluster size for chaining");
        options.addOption(REQUIRED_ANNOTATIONS, true, "Optional: string list of annotations");
        options.addOption(INDEL_ANNOTATIONS, false, "Optional: annotate clusters and TIs with INDELs");
        options.addOption(INDEL_FILE, true, "Optional: cached set of INDELs");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(LOG_VERBOSE, false, "Log extra detail");

        LinxOutput.addCmdLineArgs(options);
    }

    private static List<String> loadSampleListFile(final String filename)
    {
        try
        {
            final List<String> sampleIds = Files.readAllLines(new File(filename).toPath()).stream()
                    .filter(x -> !x.equals("SampleId"))
                    .filter(x -> !x.isEmpty())
                    .collect(Collectors.toList());

            LNX_LOGGER.info("Loaded {} specific sample IDs", sampleIds.size());

            return sampleIds;
        }
        catch (IOException exception)
        {
            LNX_LOGGER.error("failed to read sample list input CSV file({}): {}", filename, exception.toString());
            return Lists.newArrayList();
        }
    }

}
