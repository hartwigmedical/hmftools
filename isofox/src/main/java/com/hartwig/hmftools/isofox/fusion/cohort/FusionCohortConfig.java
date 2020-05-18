package com.hartwig.hmftools.isofox.fusion.cohort;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionCohortConfig
{
    private static final String FUSION_COMPARISONS = "fusion_comparisons";
    private static final String GENERATE_COHORT = "fusion_gen_cohort";
    private static final String FUSION_MIN_SAMPLES = "fusion_min_samples";
    private static final String FUSION_MIN_FRAGS = "fusion_min_frags";
    private static final String FUSION_COHORT_FILE = "fusion_cohort_file";
    private static final String WRITE_FILTERED_FUSIONS = "write_filtered_fusions";

    // parameters to include fusion in cohort file
    public final int MinSampleThreshold;
    public final int MinFragCount;
    public final boolean GenerateCohort;
    public final boolean WriteFilteredFusions;
    public final List<String> ComparisonSources;
    public final String CohortFile;

    public FusionCohortConfig(final CommandLine cmd)
    {
        MinSampleThreshold = Integer.parseInt(cmd.getOptionValue(FUSION_MIN_SAMPLES, "2"));
        MinFragCount = Integer.parseInt(cmd.getOptionValue(FUSION_MIN_FRAGS, "2"));
        GenerateCohort = cmd.hasOption(GENERATE_COHORT);
        WriteFilteredFusions = cmd.hasOption(WRITE_FILTERED_FUSIONS);

        CohortFile = cmd.getOptionValue(FUSION_COHORT_FILE);

        ComparisonSources = cmd.hasOption(FUSION_COMPARISONS) ?
                Arrays.stream(cmd.getOptionValue(FUSION_COMPARISONS).split(";")).collect(Collectors.toList()) : Lists.newArrayList();
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(GENERATE_COHORT, false, "Generate a fusion cohort file");
        options.addOption(WRITE_FILTERED_FUSIONS, false, "Apply filters to sample fusions and write to file");
        options.addOption(FUSION_MIN_SAMPLES, true, "Min number of samples to support a fusion");
        options.addOption(FUSION_MIN_FRAGS, true, "Min frag count per sample to support a fusion");
        options.addOption(FUSION_COHORT_FILE, true, "Cohort file previously generated");
        options.addOption(FUSION_COMPARISONS, true, "List of sources to compare fusions between");
    }
}
