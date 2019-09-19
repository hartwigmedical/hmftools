package com.hartwig.hmftools.bachelor.types;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.BamCountReader;
import com.hartwig.hmftools.bachelor.GermlineVcfParser;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;


public class BachelorConfig
{

    public final String RunMode;
    public final String SampleDataDir;
    public final String SampleId;
    public final String OutputDir;
    public final String PurpleDataDir;
    public final boolean IsBatchMode;
    public final int MaxBatchDirectories;
    public final List<String> RestrictedSampleIds;

    public final Map<String, Program> ProgramConfigMap;

    private boolean mIsValid;

    // config options
    public static final String CONFIG_XML = "xml_config";
    private static final String RUN_MODE = "run_mode";
    public static final String SAMPLE = "sample";
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String SAMPLE_LIST_FILE = "sample_list_file";
    private static final String BATCH_MAX_DIR = "max_batch_dir"; // only for testing

    public static final String RUN_MODE_BOTH = "Both";
    public static final String RUN_MODE_VCF_PARSE = "VcfParse";
    public static final String RUN_MODE_POST_PROCESS = "PostProcess";

    // post-process
    public static final String REF_GENOME = "ref_genome";
    private static final String PURPLE_DATA_DIRECTORY = "purple_data_dir"; // path to purple data directory
    private static final String BACH_DIRECTORY = "bachelor_dir"; // usually defaults to the 'bachelor' subdirectory of the sample dir
    private static final String BACH_INPUT_FILE = "bachelor_file"; // full path
    public static final String READ_BAMS_DIRECT = "bam_direct"; // skip BAM slicing and use of Mini-Pileup file reading

    // common
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    public static final String LOG_DEBUG = "log_debug";

    // constants
    public static final String INTERIM_FILENAME = ".bachelor_interim.csv";
    public static final String BATCH_FILE = "BATCH";

    private static final Logger LOGGER = LogManager.getLogger(BachelorConfig.class);

    public BachelorConfig(final CommandLine cmd)
    {
        mIsValid = true;

        ProgramConfigMap = Maps.newHashMap();

        if (cmd.hasOption(CONFIG_XML))
        {
            if(!loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)), ProgramConfigMap))
                mIsValid = false;
        }

        RunMode = cmd.getOptionValue(RUN_MODE, RUN_MODE_BOTH);

        LOGGER.info("Run mode: {}", RunMode);

        SampleId = cmd.getOptionValue(SAMPLE, "");

        RestrictedSampleIds = Lists.newArrayList();

        if (SampleId.isEmpty() || SampleId.equals("*"))
        {
            LOGGER.info("Running in batch mode");
            IsBatchMode = true;
            MaxBatchDirectories = Integer.parseInt(cmd.getOptionValue(BATCH_MAX_DIR, "0"));

            if (cmd.hasOption(SAMPLE_LIST_FILE))
            {
                loadSampleListFile(cmd.getOptionValue(SAMPLE_LIST_FILE), RestrictedSampleIds);
            }
        }
        else
        {
            IsBatchMode = false;
            MaxBatchDirectories = 0;
        }

        String sampleDir = cmd.getOptionValue(SAMPLE_DATA_DIR);

        if (!sampleDir.endsWith(File.separator))
        {
            sampleDir += File.separator;
        }

        SampleDataDir = sampleDir;

        String sampleOutputDir = "";

        if (cmd.hasOption(OUTPUT_DIR))
        {
            sampleOutputDir = cmd.getOptionValue(OUTPUT_DIR);

            if (!sampleOutputDir.endsWith(File.separator))
            {
                sampleOutputDir += File.separator;
            }
        }
        else
        {
            sampleOutputDir = SampleDataDir;
        }

        OutputDir = sampleOutputDir;

        PurpleDataDir = cmd.getOptionValue(PURPLE_DATA_DIRECTORY, "");
    }

    public boolean isValid() { return mIsValid; }

    public static boolean loadXML(final Path path, Map<String, Program> configMap)
    {
        try
        {
            final ConfigSchema schema = ConfigSchema.make();

            final List<Program> programs = Files.walk(path)
                    .filter(p -> p.toString().endsWith(".xml"))
                    .map(schema::processXML)
                    .filter(Objects::nonNull)
                    .collect(Collectors.toList());

            for (final Program p : programs)
            {
                if (configMap.containsKey(p.getName()))
                {
                    LOGGER.error("duplicate Programs detected: {}", p.getName());
                    return false;
                }
                else
                {
                    configMap.put(p.getName(), p);
                }
            }
        }
        catch (Exception e)
        {
            LOGGER.error("Error loading XML: {}", e.toString());
            return false;
        }

        return true;
    }

    private void loadSampleListFile(final String filename, final List<String> sampleIds)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                final String sampleId = items[0];
                sampleIds.add(sampleId);
            }

            LOGGER.info("Loaded {} specific sample ids", sampleIds.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample list input CSV file({}): {}", filename, exception.toString());
        }
    }


    @NotNull
    public static Options createOptions()
    {
        final Options options = new Options();

        // germline VCF parsing
        options.addOption(RUN_MODE, true, "VcfParse, PostProcess or Both (default)");
        options.addOption(CONFIG_XML, true, "XML with genes, black and white lists");
        options.addOption(OUTPUT_DIR, true, "When in single-sample mode, all output written to this dir");
        options.addOption(SAMPLE_DATA_DIR, true, "the run directory to look for inputs");
        options.addOption(SAMPLE_LIST_FILE, true, "Optional: limiting list of sample IDs to process");
        options.addOption(SAMPLE, true, "Sample Id (not applicable for batch mode)");
        options.addOption(BATCH_MAX_DIR, true, "Max batch directories to batch process");

        // post-process
        options.addOption(BACH_DIRECTORY, true, "Override for specific bachelor input dir, if left out then assumes in sample path & 'bachelor' sub-directory");
        options.addOption(BACH_INPUT_FILE, true, "Override for specific bachelor input file, if left out then assumes in bachelor_dir & *germline_variants.csv");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");
        options.addOption(READ_BAMS_DIRECT, false, "Read tumor alt and read depth from available BAM file");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        // logging
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        GermlineVcfParser.addCmdLineOptions(options);
        BamCountReader.addCmdLineOptions(options);

        return options;
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    public static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd)
    {
        if(!cmd.hasOption(DB_URL))
            return null;

        try
        {
            final String userName = cmd.getOptionValue(DB_USER);
            final String password = cmd.getOptionValue(DB_PASS);
            final String databaseUrl = cmd.getOptionValue(DB_URL);
            final String jdbcUrl = "jdbc:" + databaseUrl;
            return new DatabaseAccess(userName, password, jdbcUrl);
        }
        catch (SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
            return null;
        }
    }

}
