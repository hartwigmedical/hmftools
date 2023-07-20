package com.hartwig.hmftools.svtools.simulation;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ShatteringConfig
{
    public final int Iterations;
    public final int SegmentCountMin;
    public final int SegmentCountMax;
    public final boolean GroupResults;
    public final boolean RandomLinkSelection;
    public final boolean ExhuastiveSearch;

    private static final String SEG_COUNT = "sh_seg_count";
    private static final String SEG_COUNT_MIN = "sh_seg_count_min";
    private static final String SEG_COUNT_MAX = "sh_seg_count_max";
    private static final String ITERATIONS = "sh_iterations";
    private static final String COMBINE_RESULTS = "sh_group_results";
    private static final String RANDOM_SELECTION = "sh_random";
    private static final String EXHAUSTIVE_SEARCH = "sh_exhaustive";

    public ShatteringConfig(final ConfigBuilder configBuilder)
    {
        Iterations = configBuilder.getInteger(ITERATIONS);

        if(configBuilder.hasFlag(SEG_COUNT))
        {
            SegmentCountMin = configBuilder.getInteger(SEG_COUNT);
            SegmentCountMax = SegmentCountMin;
        }
        else if(configBuilder.hasFlag(SEG_COUNT_MIN) && configBuilder.hasFlag(SEG_COUNT_MAX))
        {
            SegmentCountMin = configBuilder.getInteger(SEG_COUNT_MIN);
            SegmentCountMax = configBuilder.getInteger(SEG_COUNT_MAX);
        }
        else
        {
            SegmentCountMin = SegmentCountMax = 0;
        }

        GroupResults = configBuilder.hasFlag(COMBINE_RESULTS);
        RandomLinkSelection = true; // configBuilder.hasFlag(RANDOM_SELECTION);
        ExhuastiveSearch = configBuilder.hasFlag(EXHAUSTIVE_SEARCH);
    }

    public ShatteringConfig(int segments, int iterations)
    {
        Iterations = iterations;
        SegmentCountMin = SegmentCountMax = segments;
        GroupResults = false;
        RandomLinkSelection = false;
        ExhuastiveSearch = false;
    }

    public boolean isValid()
    {
        return Iterations > 0 && SegmentCountMin> 0 && SegmentCountMax >= SegmentCountMin;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(SEG_COUNT, "Shattering segments", 1);
        configBuilder.addInteger(SEG_COUNT_MIN, "Shattering segments range min", 1);
        configBuilder.addInteger(SEG_COUNT_MAX, "Shattering segments range max", 1);
        configBuilder.addInteger(ITERATIONS, "Shattering test iterations", 1);
        configBuilder.addFlag(COMBINE_RESULTS, "Shattering group like results");
        configBuilder.addFlag(RANDOM_SELECTION, "Shattering use random selection of next link");
        configBuilder.addFlag(EXHAUSTIVE_SEARCH, "Shattering find all possible link combinations");
    }

}
