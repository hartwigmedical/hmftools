package com.hartwig.hmftools.purple.config;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class ChartConfig
{
    private static final String CIRCOS = "circos";
    private static final String DISABLE = "no_charts";

    public final boolean Disabled;

    public final String CircosBinary;

    public final String PlotDirectory;
    public final String CircosDirectory;

    public ChartConfig(final CommandLine cmd, final String outputDir)
    {
        Disabled = cmd.hasOption(DISABLE);
        PlotDirectory = outputDir + File.separator + "plot";
        CircosDirectory = outputDir + File.separator + "circos";
        CircosBinary = cmd.getOptionValue(CIRCOS);
    }

    public static void addOptions(@NotNull Options options)
    {
        options.addOption(CIRCOS, true, "Location of circos binary");
        options.addOption(DISABLE, false, "Disable charts");
    }

}
