package com.hartwig.hmftools.isofox.fusion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionConfig
{
    public final boolean WriteChimericReads;
    public final String ReadsFile;
    public final boolean PerformanceStats;
    public final boolean RecoverMissingReads;

    private static final String WRITE_CHIMERIC_READS = "write_chimeric_reads";
    private static final String CHIMERIC_READS_FILE = "chimeric_reads_file";
    private static final String RECOVER_MISSING_READS = "recover_chimeric_reads";

    public FusionConfig(final CommandLine cmd)
    {
        WriteChimericReads = cmd.hasOption(WRITE_CHIMERIC_READS);
        ReadsFile = cmd.getOptionValue(CHIMERIC_READS_FILE);
        RecoverMissingReads = false;
        PerformanceStats = true;
    }

    public FusionConfig()
    {
        WriteChimericReads = false;
        ReadsFile = null;
        RecoverMissingReads = false;
        PerformanceStats = false;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(RECOVER_MISSING_READS, false, "Research BAM for incomplete supplementary reads");
        options.addOption(WRITE_CHIMERIC_READS, false, "Write chimeric read data");
        options.addOption(CHIMERIC_READS_FILE, true, "Chimeric reads to use instead of extracting from BAM");
    }

}
