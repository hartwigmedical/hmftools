package com.hartwig.hmftools.purple;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
        PlotDirectory = outputDir + "plot" + File.separator;
        CircosDirectory = outputDir + "circos" + File.separator;
        CircosBinary = configBuilder.getValue(CIRCOS);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CIRCOS, false, "Location of circos binary");
        configBuilder.addFlag(DISABLE, "Disable charts");
    }

}
