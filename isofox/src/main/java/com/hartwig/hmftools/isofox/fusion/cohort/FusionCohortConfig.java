package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.isofox.fusion.FusionConfig.FUSION_COHORT_FILE;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FusionCohortConfig
{
    private static final String FUSION_COMPARISONS = "fusion_comparisons";
    private static final String GENERATE_COHORT = "fusion_gen_cohort";
    private static final String FUSION_MIN_SAMPLES = "fusion_min_samples";
    private static final String FUSION_MIN_FRAGS = "fusion_min_frags";
    private static final String WRITE_FILTERED_FUSIONS = "fusion_write_filtered";
    private static final String WRITE_COMBINED_FUSIONS = "fusion_write_combined";
    private static final String COMPARE_UNFILTERED = "fusion_compare_unfiltered";
    private static final String REWRITE_ANNOTATED = "fusion_rewrite_annotated";

    // cohort file generation and filters
    public final boolean GenerateCohort;
    public final int MinSampleThreshold;
    public final int MinFragCount;

    public final String CohortFile;

    public final boolean WriteFilteredFusions; // rewrite passing fusions
    public final boolean WriteCombinedFusions;
    public final boolean RewriteAnnotatedFusions; // rewrite all fusions with cohort and other gene data

    public final String ComparisonSource;
    public final boolean CompareUnfiltered;

    public FusionCohortConfig(final CommandLine cmd)
    {
        MinSampleThreshold = Integer.parseInt(cmd.getOptionValue(FUSION_MIN_SAMPLES, "2"));
        MinFragCount = Integer.parseInt(cmd.getOptionValue(FUSION_MIN_FRAGS, "2"));
        GenerateCohort = cmd.hasOption(GENERATE_COHORT);
        WriteFilteredFusions = cmd.hasOption(WRITE_FILTERED_FUSIONS);
        CompareUnfiltered = cmd.hasOption(COMPARE_UNFILTERED);
        RewriteAnnotatedFusions = cmd.hasOption(REWRITE_ANNOTATED);
        WriteCombinedFusions = cmd.hasOption(WRITE_COMBINED_FUSIONS);

        CohortFile = cmd.getOptionValue(FUSION_COHORT_FILE);
        ComparisonSource = cmd.getOptionValue(FUSION_COMPARISONS);
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(GENERATE_COHORT, false, "Generate a fusion cohort file");
        options.addOption(WRITE_FILTERED_FUSIONS, false, "Apply filters to sample fusions and write to file");
        options.addOption(FUSION_MIN_SAMPLES, true, "Min number of samples to support a fusion");
        options.addOption(FUSION_MIN_FRAGS, true, "Min frag count per sample to support a fusion");
        options.addOption(FUSION_COHORT_FILE, true, "Cohort file previously generated");
        options.addOption(FUSION_COMPARISONS, true, "List of sources to compare fusions between");
        options.addOption(COMPARE_UNFILTERED, false, "Included unfiltered fusions in comparison with external tools");
        options.addOption(REWRITE_ANNOTATED, false, "Rewrtite all fusions with cohort and other annotations");
        options.addOption(WRITE_COMBINED_FUSIONS, false, "Write a combined file with passing fusion");

        options.addOption(KNOWN_FUSIONS_FILE, true, "Known fusion file");
    }
}
