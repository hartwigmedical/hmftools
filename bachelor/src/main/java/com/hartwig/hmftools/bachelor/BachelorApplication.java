package com.hartwig.hmftools.bachelor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.types.ConfigSchema;

import nl.hartwigmedicalfoundation.bachelor.Program;

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
    ExternalDBFilters mFilterFileBuilder;
    private Map<String, Program> mConfigMap;

    private static final String RUN_MODE_BOTH = "Both";
    private static final String RUN_MODE_VCF_PARSE = "VcfParse";
    private static final String RUN_MODE_POST_PROCESS = "PostProcess";

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);

    private BachelorApplication()
    {
        mGermlineVcfParser = null;
        mPostProcessor = null;
        mFilterFileBuilder = null;
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        if (cmd.hasOption(CONFIG_XML))
        {
            try
            {
                mConfigMap = loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
            }
            catch(Exception e)
            {
                LOGGER.error("error loading XML: {}", e.toString());
                return false;
            }
        }

        String outputDir = cmd.getOptionValue(OUTPUT_DIR);

        if(cmd.hasOption(CREATE_FILTER_FILE))
        {
            LOGGER.info("building Clinvar filter files");
            final String filterInputFile = cmd.getOptionValue(CREATE_FILTER_FILE);
            mFilterFileBuilder = new ExternalDBFilters(filterInputFile, outputDir);
            return true;
        }

        String sample = cmd.getOptionValue(SAMPLE);
        List<String> sampleIds = Lists.newArrayList();

        if (sample == null || sample.equals("*"))
        {
            LOGGER.info("running in batch mode");

            if(cmd.hasOption(SAMPLE_LIST_FILE))
            {
                sampleIds = loadSampleListFile(cmd.getOptionValue(SAMPLE_LIST_FILE));
            }
        }
        else
        {
            sampleIds.add(sample);
        }

        String sampleDataDirectory = cmd.getOptionValue(SAMPLE_DATA_DIR);

        if (!sampleDataDirectory.endsWith(File.separator))
        {
            sampleDataDirectory += File.separator;
        }

        String runMode = cmd.getOptionValue(RUN_MODE, RUN_MODE_BOTH);

        if(runMode == RUN_MODE_BOTH || runMode == RUN_MODE_VCF_PARSE)
        {
            mGermlineVcfParser = new GermlineVcfParser();
            mGermlineVcfParser.initialise(cmd, mConfigMap, sampleIds, sampleDataDirectory, outputDir);
        }

        if(runMode == RUN_MODE_BOTH || runMode == RUN_MODE_POST_PROCESS)
        {
            mPostProcessor = new BachelorPostProcess();
            mPostProcessor.initialise(cmd, sampleIds, sampleDataDirectory);
        }

        return true;

    }

    public void run()
    {
        if(mFilterFileBuilder != null)
        {
            final Program program = mConfigMap.values().iterator().next();

            mFilterFileBuilder.createFilterFile(program);
            LOGGER.info("run complete");
        }

        if(mGermlineVcfParser != null)
        {
            mGermlineVcfParser.run();
        }

        if(mPostProcessor != null)
        {
            if(mGermlineVcfParser != null)
                mPostProcessor.setBachelorRecords(mGermlineVcfParser.getBachelorRecords());
            else
                mPostProcessor.loadBachelorRecords();

            mPostProcessor.run();
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

    private static Map<String, Program> loadXML(final Path path) throws IOException, SAXException
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

    private static final String CONFIG_XML = "configXml";
    private static final String RUN_MODE = "run_mode";
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String SAMPLE = "sample";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SAMPLE_LIST_FILE = "sample_list_file";
    private static final String CREATE_FILTER_FILE = "create_filter_file";

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(CONFIG_XML, true, "single config XML to run");
        options.addOption(OUTPUT_DIR, true, "output file");
        options.addOption(SAMPLE_DATA_DIR, true, "the run directory to look for inputs");
        options.addOption(SAMPLE_LIST_FILE, true, "Optional: limiting list of sample IDs to process");
        options.addOption(SAMPLE, true, "sample id");
        options.addOption(CREATE_FILTER_FILE, true, "Optional: create black and white list filter files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        GermlineVcfParser.addCmdLineOptions(options);
        BachelorPostProcess.addCmdLineOptions(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException
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
