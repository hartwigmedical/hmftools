package com.hartwig.hmftools.svanalysis.simulation;

import static com.hartwig.hmftools.svanalysis.SvAnalyser.DATA_OUTPUT_PATH;

import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SvSimulator
{
    private SvSimShattering mSimShattering;

    public SvSimulator()
    {
        mSimShattering = new SvSimShattering();
    }

    public static void addCmdLineArgs(Options options)
    {

    }

    public void loadConfig(final CommandLine cmd)
    {
        mSimShattering.setOutputDir(cmd.getOptionValue(DATA_OUTPUT_PATH));
    }

    public void run()
    {
        mSimShattering.initialise(4, 1);
        mSimShattering.run();
        mSimShattering.logResults();
        mSimShattering.close();
    }

}

