package com.hartwig.hmftools.linx.simulation;

import com.hartwig.hmftools.linx.types.SvaConfig;

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

public class SvSimulator
{
    private SvSimShattering mSimShattering;

    private static String SHATTERING_SEG_COUNT = "sim_sh_seg_count";
    private static String SHATTERING_ITERATIONS = "sim_sh_iterations";

    private static final Logger LOGGER = LogManager.getLogger(SvSimulator.class);

    public SvSimulator()
    {
        mSimShattering = new SvSimShattering();
    }

    public void loadConfig(final CommandLine cmd, final String outputDir)
    {
        mSimShattering.setOutputDir(outputDir);

        int segCount = Integer.parseInt(cmd.getOptionValue(SHATTERING_SEG_COUNT, "1"));
        int iterations = Integer.parseInt(cmd.getOptionValue(SHATTERING_ITERATIONS, "1"));
        mSimShattering.initialise(segCount, iterations);
    }

    public void run()
    {
        LOGGER.info("starting simulations");
        mSimShattering.run();
        mSimShattering.logResults();
        mSimShattering.close();
        LOGGER.info("simulations complete");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SHATTERING_SEG_COUNT, true, "Shattering segments");
        options.addOption(SHATTERING_ITERATIONS, true, "Shattering iterations");
        options.addOption(SvaConfig.DATA_OUTPUT_DIR, true, "Output directory");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(SvaConfig.LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String outputDir = SvaConfig.formOutputPath(cmd.getOptionValue(SvaConfig.DATA_OUTPUT_DIR));
        SvSimulator simulator = new SvSimulator();
        simulator.loadConfig(cmd, outputDir);
        simulator.run();
    }

}

