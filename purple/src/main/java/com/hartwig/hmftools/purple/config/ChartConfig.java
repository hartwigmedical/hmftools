package com.hartwig.hmftools.purple.config;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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

    public ChartConfig(final ConfigBuilder configBuilder, final String outputDir)
    {
        Disabled = configBuilder.hasFlag(DISABLE);
        PlotDirectory = outputDir + "plot";
        CircosDirectory = outputDir + "circos";
        CircosBinary = configBuilder.getValue(CIRCOS);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPathItem(CIRCOS, false, "Location of circos binary");
        configBuilder.addFlagItem(DISABLE, "Disable charts");
    }

}
