package com.hartwig.hmftools.bachelor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.bachelor.types.ConfigSchema;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.xml.sax.SAXException;

public class BachelorApplication {

    GermlineVcfParser mGermlineVcfParser;
    BachelorPostProcess mPostProcessor;

    private Map<String, Program> mConfigMap;
    private String mSampleDataDir;
    private String mSingleSampleId;
    private List<String> mRestrictedSampleIds;
    private List<RunDirectory> mSampleDataDirectories;
    private boolean mIsBatchMode;
    private int mMaxBatchDirectories;

    // config options
    public static final String CONFIG_XML = "xml_config";
    public static final String BATCH_OUTPUT_DIR = "batch_output_dir";
    public static final String LOG_DEBUG = "log_debug";

    private static final String RUN_MODE = "run_mode";
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String SAMPLE = "sample";
    private static final String SAMPLE_LIST_FILE = "sample_list_file";
    private static final String BATCH_MAX_DIR = "max_batch_dir"; // only for testing

    private static final String RUN_MODE_BOTH = "Both";
    private static final String RUN_MODE_VCF_PARSE = "VcfParse";
    private static final String RUN_MODE_POST_PROCESS = "PostProcess";

    public static final String DEFAULT_BACH_DIRECTORY = "bachelor";

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);

    private BachelorApplication()
    {
        mGermlineVcfParser = null;
        mPostProcessor = null;
        mSampleDataDir = "";
        mSingleSampleId = "";
        mSampleDataDirectories = Lists.newArrayList();
        mRestrictedSampleIds = Lists.newArrayList();
        mIsBatchMode = false;
        mMaxBatchDirectories = 0;
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        if (!cmd.hasOption(CONFIG_XML))
            return false;

        try
        {
            mConfigMap = loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
        }
        catch(Exception e)
        {
            LOGGER.error("error loading XML: {}", e.toString());
            return false;
        }

        String runMode = cmd.getOptionValue(RUN_MODE, RUN_MODE_BOTH);

        LOGGER.info("run mode: {}", runMode);

        String batchOutputDir = cmd.getOptionValue(BATCH_OUTPUT_DIR, "");

        mSingleSampleId = cmd.getOptionValue(SAMPLE, "");

        if (mSingleSampleId.isEmpty() || mSingleSampleId.equals("*"))
        {
            LOGGER.info("running in batch mode");
            mIsBatchMode = true;

            if(cmd.hasOption(SAMPLE_LIST_FILE))
            {
                mRestrictedSampleIds = loadSampleListFile(cmd.getOptionValue(SAMPLE_LIST_FILE));
            }
        }

        mMaxBatchDirectories = Integer.parseInt(cmd.getOptionValue(BATCH_MAX_DIR, "0"));

        mSampleDataDir = cmd.getOptionValue(SAMPLE_DATA_DIR);

        if (!mSampleDataDir.endsWith(File.separator))
        {
            mSampleDataDir += File.separator;
        }

        setSampleDataDirectories();

        if(runMode.equals(RUN_MODE_BOTH) || runMode.equals(RUN_MODE_VCF_PARSE))
        {
            mGermlineVcfParser = new GermlineVcfParser();
            mGermlineVcfParser.initialise(cmd, mConfigMap, mIsBatchMode, batchOutputDir);
        }

        if(runMode.equals(RUN_MODE_BOTH) || runMode.equals(RUN_MODE_POST_PROCESS))
        {
            mPostProcessor = new BachelorPostProcess();
            mPostProcessor.initialise(cmd, mIsBatchMode, batchOutputDir);
        }

        return true;
    }

    public void run()
    {
        for (int i = 0; i < mSampleDataDirectories.size(); ++i)
        {
            final RunDirectory runDir = mSampleDataDirectories.get(i);

            String sampleId = "";

            if(!mIsBatchMode)
            {
                sampleId = mSingleSampleId;
            }
            else
            {
                try
                {
                    final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDir.sampleDir().toString());
                    sampleId = runContext.tumorSample();
                }
                catch (Exception e)
                {
                    // Skip using meta data
                    // sampleId = runDir.getPatientID();
                    continue;
                }
            }

            if(!mRestrictedSampleIds.isEmpty() && !mRestrictedSampleIds.contains(sampleId))
            {
                LOGGER.info("skipping sampleId({}) not in specified list", sampleId);
                continue;
            }

            // processSampleDirectory(runDir, "");

            if(mGermlineVcfParser != null)
            {
                mGermlineVcfParser.run(runDir, sampleId);
            }

            if(mPostProcessor != null)
            {
                mPostProcessor.run(runDir, sampleId, mGermlineVcfParser != null ? mGermlineVcfParser.getBachelorRecords() : null);
            }

            if(mMaxBatchDirectories > 0 && i >= mMaxBatchDirectories)
                break;
        }

        if(mGermlineVcfParser != null)
        {
            mGermlineVcfParser.close();
        }

        if(mPostProcessor != null)
        {
            mPostProcessor.close();
        }

        LOGGER.info("run complete");
    }

    private void setSampleDataDirectories()
    {
        final Path sampleDataPath = Paths.get(mSampleDataDir);

        if(mIsBatchMode)
        {
            try (final Stream<Path> stream = Files.walk(sampleDataPath, 1, FileVisitOption.FOLLOW_LINKS).parallel())
            {
                mSampleDataDirectories = stream.filter(p -> p.toFile().isDirectory())
                        .filter(p -> !p.equals(sampleDataPath))
                        .map(RunDirectory::new)
                        .collect(Collectors.toList());
            }
            catch (Exception e)
            {
                LOGGER.error("failed walking batch directories: {}", e.toString());
            }

            LOGGER.info("found {} batch directories", mSampleDataDirectories.size());
        }
        else
        {
            mSampleDataDirectories.add(new RunDirectory(sampleDataPath));
        }
    }

    private List<String> loadSampleListFile(final String filename)
    {
        List<String> sampleIds = Lists.newArrayList();

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

        return sampleIds;
    }

    public static Map<String, Program> loadXML(final Path path) throws IOException, SAXException
    {
        final ConfigSchema schema = ConfigSchema.make();

        final List<Program> programs = Files.walk(path)
                .filter(p -> p.toString().endsWith(".xml"))
                .map(schema::processXML)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        final Map<String, Program> result = Maps.newHashMap();

        for (final Program p : programs)
        {
            if (result.containsKey(p.getName()))
            {
                LOGGER.error("duplicate Programs detected: {}", p.getName());
                System.exit(1);
            }
            else
            {
                result.put(p.getName(), p);
            }
        }

        return result;
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(RUN_MODE, true, "VcfParse, PostProcess or Both (default)");
        options.addOption(CONFIG_XML, true, "XML with genes, black and white lists");
        options.addOption(BATCH_OUTPUT_DIR, true, "Optional: when in batch mode, all output written to single file");
        options.addOption(SAMPLE_DATA_DIR, true, "the run directory to look for inputs");
        options.addOption(SAMPLE_LIST_FILE, true, "Optional: limiting list of sample IDs to process");
        options.addOption(SAMPLE, true, "Sample Id (not applicable for batch mode)");
        options.addOption(BATCH_MAX_DIR, true, "Max batch directories to batch process");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        GermlineVcfParser.addCmdLineOptions(options);
        BachelorPostProcess.addCmdLineOptions(options);
        BamCountReader.addCmdLineOptions(options);

        return options;
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public static void main(final String... args)
    {
        final Options options = createOptions();

        BachelorApplication bachelorApp = new BachelorApplication();

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            if(!bachelorApp.loadConfig(cmd))
            {
                System.exit(1);
                return;
            }

            bachelorApp.run();
        }
        catch (final ParseException e)
        {
            printHelpAndExit(options);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    private static void printHelpAndExit(final Options options)
    {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bachelor", "Determines eligibility", options, "", true);
        System.exit(1);
    }

}
