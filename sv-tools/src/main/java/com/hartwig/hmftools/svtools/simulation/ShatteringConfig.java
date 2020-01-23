package com.hartwig.hmftools.svtools.simulation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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

    public ShatteringConfig(final CommandLine cmd)
    {
        Iterations = Integer.parseInt(cmd.getOptionValue(ITERATIONS, "1"));

        if(cmd.hasOption(SEG_COUNT))
        {
            SegmentCountMin = Integer.parseInt(cmd.getOptionValue(SEG_COUNT, "1"));
            SegmentCountMax = SegmentCountMin;
        }
        else if(cmd.hasOption(SEG_COUNT_MIN) && cmd.hasOption(SEG_COUNT_MAX))
        {
            SegmentCountMin = Integer.parseInt(cmd.getOptionValue(SEG_COUNT_MIN, "1"));
            SegmentCountMax = Integer.parseInt(cmd.getOptionValue(SEG_COUNT_MAX, "1"));
        }
        else
        {
            SegmentCountMin = SegmentCountMax = 0;
        }

        GroupResults = cmd.hasOption(COMBINE_RESULTS);
        RandomLinkSelection = true; // cmd.hasOption(RANDOM_SELECTION);
        ExhuastiveSearch = cmd.hasOption(EXHAUSTIVE_SEARCH);
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

    public static void addCommandLineOptions(Options options)
    {
        options.addOption(SEG_COUNT, true, "Shattering segments");
        options.addOption(SEG_COUNT_MIN, true, "Shattering segments range min");
        options.addOption(SEG_COUNT_MAX, true, "Shattering segments range max");
        options.addOption(ITERATIONS, true, "Shattering test iterations");
        options.addOption(COMBINE_RESULTS, false, "Shattering group like results");
        options.addOption(RANDOM_SELECTION, false, "Shattering use random selection of next link");
        options.addOption(EXHAUSTIVE_SEARCH, false, "Shattering find all possible link combinations");
    }

}
