package com.hartwig.hmftools.sig_analyser;

import com.hartwig.hmftools.sig_analyser.buckets.BucketAnalyser;
import com.hartwig.hmftools.sig_analyser.nmf.NmfConfig;
import com.hartwig.hmftools.sig_analyser.nmf.NmfManager;
import com.hartwig.hmftools.sig_analyser.sim.SampleSimulator;
import com.hartwig.hmftools.sig_analyser.sim.SimConfig;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SigAnalyser
{
    public static final String GENERIC_INPUT_FILE = "gen_input_file";
    public static final String LOG_DEBUG = "log_debug";

    public static final String OUTPUT_DIR = "output_dir";
    public static final String OUTPUT_FILE_ID = "output_file_id";

    private static final String RUN_NMF = "run_nmf";
    private static final String RUN_BA = "run_buckets";
    private static final String RUN_SAMPLE_SIM = "run_sim";

    public static final Logger SIG_LOGGER = LogManager.getLogger(SigAnalyser.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = createBasicOptions();

        // allow components to add their own arg lists
        SimConfig.addCmdLineArgs(options);
        NmfConfig.addCmdLineArgs(options);
        BucketAnalyser.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        SIG_LOGGER.info("starting signature analyser");

        final GenericDataCollection collection = GenericDataLoader.loadFile(cmd.getOptionValue(GENERIC_INPUT_FILE));

        if(cmd.hasOption(RUN_NMF))
        {

            NmfManager nmfManager = new NmfManager();
            nmfManager.initialise(collection, cmd);
            nmfManager.run();
            // nmfManager.runTests();
        }

        if(cmd.hasOption(RUN_BA))
        {
            BucketAnalyser bucketAnalyser = new BucketAnalyser(collection, cmd);
            bucketAnalyser.run();
        }

        if(cmd.hasOption(RUN_SAMPLE_SIM))
        {
            SampleSimulator sampleSimulator = new SampleSimulator();
            sampleSimulator.initialise(cmd);
            sampleSimulator.run();
        }

        SIG_LOGGER.info("analysis complete");
    }

    @NotNull
    private static Options createBasicOptions()
    {
        Options options = new Options();
        options.addOption(GENERIC_INPUT_FILE, true, "Path to the main input file");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");

        options.addOption(RUN_NMF, false, "Run NMF");
        options.addOption(RUN_BA, false, "Run bucket analysis");
        options.addOption(RUN_SAMPLE_SIM, false, "Generate simulated sample counts");

        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}