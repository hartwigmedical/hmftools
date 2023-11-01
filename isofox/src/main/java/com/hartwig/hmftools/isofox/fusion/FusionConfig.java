package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.DEFAULT_HARD_FILTER_MIN_FRAGS;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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

    public FusionConfig(final ConfigBuilder configBuilder)
    {
        WriteChimericReads = configBuilder.hasFlag(WRITE_CHIMERIC_READS);
        WriteChimericFragments = configBuilder.hasFlag(WRITE_CHIMERIC_FRAGS);
        RunPerfChecks = configBuilder.hasFlag(RUN_FUSION_PERF);
        CohortFile = configBuilder.getValue(FUSION_COHORT_FILE);
        MinHardFilterFrags = configBuilder.getInteger(MIN_FRAGS_HARD_FILTER);

        KnownFusions = new KnownFusionCache();
        KnownFusions.loadFromFile(configBuilder);

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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(WRITE_CHIMERIC_READS, "Write chimeric read data");
        configBuilder.addFlag(WRITE_CHIMERIC_FRAGS, "Write chimeric fragment data");
        addKnownFusionFileOption(configBuilder);
        configBuilder.addPath(FUSION_COHORT_FILE, false, "Cohort file previously generated");
        configBuilder.addInteger(MIN_FRAGS_HARD_FILTER, "Hard filter chimeric translocations", DEFAULT_HARD_FILTER_MIN_FRAGS);
        configBuilder.addFlag(RUN_FUSION_PERF, "Write chimeric fragment data");
    }
}
