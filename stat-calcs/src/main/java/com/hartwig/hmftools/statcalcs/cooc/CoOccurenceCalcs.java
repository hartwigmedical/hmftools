package com.hartwig.hmftools.statcalcs.cooc;

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

public class CoOccurenceCalcs
{
    private ThreeVarCoOccurence mThreeVarCoOccurence;
    private TwoVarCoOccurence mTwoVarCoOccurence;
    private SampleCountsCoOccurence mSampleCountsCoOccurence;

    private static final String LOG_DEBUG = "log_debug";
    private static final String DATA_OUTPUT_DIR = "output_dir";

    private static final Logger LOGGER = LogManager.getLogger(CoOccurenceCalcs.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String outputDir = cmd.getOptionValue(DATA_OUTPUT_DIR);

        CoOccurenceCalcs statsRoutines = new CoOccurenceCalcs();
        statsRoutines.loadConfig(cmd, outputDir);
        statsRoutines.runStatistics();
        LOGGER.info("run complete");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Verbose logging");

        TwoVarCoOccurence.addCmdLineOptions(options);
        ThreeVarCoOccurence.addCmdLineOptions(options);
        SampleCountsCoOccurence.addCmdLineOptions(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public CoOccurenceCalcs()
    {
        mThreeVarCoOccurence = null;
        mSampleCountsCoOccurence = null;
        mTwoVarCoOccurence = null;
    }

    public boolean loadConfig(final CommandLine cmd, final String outputDir)
    {
        boolean valid = true;

        if(SampleCountsCoOccurence.hasConfig(cmd))
        {
            mSampleCountsCoOccurence = new SampleCountsCoOccurence(cmd, outputDir);
        }

        if(ThreeVarCoOccurence.hasConfig(cmd))
        {
            mThreeVarCoOccurence = new ThreeVarCoOccurence(cmd, outputDir);
        }

        if(TwoVarCoOccurence.hasConfig(cmd))
        {
            mTwoVarCoOccurence = new TwoVarCoOccurence(cmd, outputDir);
        }

        return valid;
    }

    public void runStatistics()
    {
        if(mSampleCountsCoOccurence != null)
            mSampleCountsCoOccurence.run();

        if(mThreeVarCoOccurence != null)
            mThreeVarCoOccurence.run();

        if(mTwoVarCoOccurence != null)
            mTwoVarCoOccurence.run();
    }

}


