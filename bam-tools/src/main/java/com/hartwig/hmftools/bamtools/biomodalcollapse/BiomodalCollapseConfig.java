package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BiomodalCollapseConfig
{
    public final String Fastq1Path;
    public final String Fastq2Path;
    public final String CollapsedFastqOutputPath;
    public final int Threads;
    public final int MaxFastqPairsProcessed;
    public final String RefResolvedFastqPath;
    public final String DebugStatsOutputPath;

    private static final String FASTQ1_PATH = "fastq1";
    private static final String FASTQ2_PATH = "fastq2";
    private static final String COLLAPSED_FASTQ_OUTPUT_PATH = "out_fastq";
    private static final String MAX_FASTQ_PAIRS_PROCESSED = "max_fastq_pairs_read";
    private static final String REF_RESOLVED_FASTQ_PATH = "ref_resolved_fastq";
    private static final String DEBUG_STATS_OUTPUT_PATH = "out_debug_tsv";

    public BiomodalCollapseConfig(final ConfigBuilder configBuilder)
    {

        Fastq1Path = configBuilder.getValue(FASTQ1_PATH);
        Fastq2Path = configBuilder.getValue(FASTQ2_PATH);
        CollapsedFastqOutputPath = configBuilder.getValue(COLLAPSED_FASTQ_OUTPUT_PATH);
        Threads = max(parseThreads(configBuilder), 1);
        MaxFastqPairsProcessed = configBuilder.getInteger(MAX_FASTQ_PAIRS_PROCESSED);

        RefResolvedFastqPath = configBuilder.getValue(REF_RESOLVED_FASTQ_PATH, null);
        DebugStatsOutputPath = configBuilder.getValue(DEBUG_STATS_OUTPUT_PATH, null);

        validate();
    }

    private void validate()
    {
        if(RefResolvedFastqPath != null && DebugStatsOutputPath == null)
        {
            throw new RuntimeException(format("Cannot specify %s unless %s is specified", REF_RESOLVED_FASTQ_PATH, DEBUG_STATS_OUTPUT_PATH));
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(FASTQ1_PATH, true, "Path to the fastq file of first in pair reads");
        configBuilder.addConfigItem(FASTQ2_PATH, true, "Path to the fastq file of second in pair reads");
        configBuilder.addConfigItem(COLLAPSED_FASTQ_OUTPUT_PATH, true, "Path where the final collapsed fastq is written to");
        configBuilder.addInteger(MAX_FASTQ_PAIRS_PROCESSED, "Once this many fastq pairs are processed we exit", -1);

        configBuilder.addConfigItem(REF_RESOLVED_FASTQ_PATH, false, "Path to resolved fragments for comparison");
        configBuilder.addConfigItem(DEBUG_STATS_OUTPUT_PATH, false, "Path to where debug details for each read pair will be written to as a tsv");

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
