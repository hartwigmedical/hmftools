package com.hartwig.hmftools.svanalysis.simulation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SHATTERING_SEG_COUNT, true, "Sim: shattering segments");
        options.addOption(SHATTERING_ITERATIONS, true, "Sim: shattering iterations");
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

}

