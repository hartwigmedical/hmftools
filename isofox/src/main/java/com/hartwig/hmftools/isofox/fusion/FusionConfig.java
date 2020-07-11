package com.hartwig.hmftools.isofox.fusion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionConfig
{
    public final boolean WriteChimericReads;
    public final boolean PerformanceStats;

    private static final String WRITE_CHIMERIC_READS = "write_chimeric_reads";

    public FusionConfig(final CommandLine cmd)
    {
        WriteChimericReads = cmd.hasOption(WRITE_CHIMERIC_READS);
        PerformanceStats = true;
    }

    public FusionConfig()
    {
        WriteChimericReads = false;
        PerformanceStats = false;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(WRITE_CHIMERIC_READS, false, "Write chimeric read data");
    }

}
