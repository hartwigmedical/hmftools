package com.hartwig.hmftools.linx.simulation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class ShatteringConfig
{

    public final int Iterations;
    public final int SegmentCountMin;
    public final int SegmentCountMax;

    private static String SEG_COUNT = "shattering_seg_count";
    private static String SEG_COUNT_MIN = "shattering_seg_count_min";
    private static String SEG_COUNT_MAX = "shattering_seg_count_max";
    private static String ITERATIONS = "shattering_iterations";
    // private static String ITERATIONS = "shattering_iterations";

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
    }

    public ShatteringConfig(int segments, int iterations)
    {
        Iterations = iterations;
        SegmentCountMin = SegmentCountMax = segments;
    }

    public boolean isValid()
    {
        return Iterations > 0 && SegmentCountMin> 0 && SegmentCountMax >= SegmentCountMin;
    }

    public static void addCommandLineOptions(Options options)
    {
        options.addOption(SEG_COUNT, true, "Shattering segments");
        options.addOption(SEG_COUNT_MIN, true, "Shattering segments");
        options.addOption(SEG_COUNT_MAX, true, "Shattering segments");
        options.addOption(ITERATIONS, true, "Shattering iterations");
    }

}
