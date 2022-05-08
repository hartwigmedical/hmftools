package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.DEFAULT_HARD_FILTER_MIN_FRAGS;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionConfig
{
    public final boolean WriteChimericReads;
    public final boolean WriteChimericFragments;
    public final boolean PerformanceStats;
    public final String CohortFile;
    public final boolean CacheFragments;
    public final boolean RunPerfChecks;
    public int MinHardFilterFrags;

    public final KnownFusionCache KnownFusions;

    private static final String WRITE_CHIMERIC_READS = "write_chimeric_reads";
    private static final String WRITE_CHIMERIC_FRAGS = "write_chimeric_frags";
    private static final String MIN_FRAGS_HARD_FILTER = "fusion_min_frags_filter";
    private static final String RUN_FUSION_PERF = "run_fusion_perfs";

    public static final String FUSION_COHORT_FILE = "fusion_cohort_file";

    public FusionConfig(final CommandLine cmd)
    {
        WriteChimericReads = cmd.hasOption(WRITE_CHIMERIC_READS);
        WriteChimericFragments = cmd.hasOption(WRITE_CHIMERIC_FRAGS);
        RunPerfChecks = cmd.hasOption(RUN_FUSION_PERF);
        CohortFile = cmd.getOptionValue(FUSION_COHORT_FILE);
        MinHardFilterFrags = Integer.parseInt(cmd.getOptionValue(MIN_FRAGS_HARD_FILTER, String.valueOf(DEFAULT_HARD_FILTER_MIN_FRAGS)));

        KnownFusions = new KnownFusionCache();
        KnownFusions.loadFromFile(cmd);

        CacheFragments = WriteChimericFragments || WriteChimericReads;

        PerformanceStats = true;
    }

    public FusionConfig()
    {
        WriteChimericReads = false;
        WriteChimericFragments = false;
        CacheFragments = true;
        KnownFusions = new KnownFusionCache();
        PerformanceStats = false;
        CohortFile = null;
        RunPerfChecks = false;
        MinHardFilterFrags = 0;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(WRITE_CHIMERIC_READS, false, "Write chimeric read data");
        options.addOption(WRITE_CHIMERIC_FRAGS, false, "Write chimeric fragment data");
        options.addOption(KNOWN_FUSIONS_FILE, true, KNOWN_FUSIONS_FILE_DESC);
        options.addOption(FUSION_COHORT_FILE, true, "Cohort file previously generated");
        options.addOption(MIN_FRAGS_HARD_FILTER, true, "Hard filter chimeric translocations");
        options.addOption(RUN_FUSION_PERF, false, "Write chimeric fragment data");
    }
}
