package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionConfig
{
    public final boolean WriteChimericReads;
    public final boolean WriteChimericFragments;
    public final boolean PerformanceStats;
    public final String ChimericReadsFile;

    public final KnownFusionCache KnownFusions;

    private static final String WRITE_CHIMERIC_READS = "write_chimeric_reads";
    private static final String WRITE_CHIMERIC_FRAGS = "write_chimeric_frags";
    private static final String CHIMERIC_READ_FILE = "chimeric_reads_file";

    public FusionConfig(final CommandLine cmd)
    {
        WriteChimericReads = cmd.hasOption(WRITE_CHIMERIC_READS);
        WriteChimericFragments = cmd.hasOption(WRITE_CHIMERIC_FRAGS);
        ChimericReadsFile = cmd.getOptionValue(CHIMERIC_READ_FILE);
        KnownFusions = new KnownFusionCache();
        KnownFusions.loadFromFile(cmd);

        PerformanceStats = true;
    }

    public FusionConfig()
    {
        WriteChimericReads = false;
        WriteChimericFragments = false;
        ChimericReadsFile = null;
        PerformanceStats = false;
        KnownFusions = new KnownFusionCache();
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(WRITE_CHIMERIC_READS, false, "Write chimeric read data");
        options.addOption(WRITE_CHIMERIC_FRAGS, false, "Write chimeric fragment data");
        options.addOption(CHIMERIC_READ_FILE, true, "Replay chimeric reads from file without BAM");
        options.addOption(KNOWN_FUSIONS_FILE, true, "Known fusion file");
    }
}
