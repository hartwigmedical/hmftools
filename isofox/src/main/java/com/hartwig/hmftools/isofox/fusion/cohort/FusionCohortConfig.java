package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.isofox.fusion.FusionConfig.FUSION_COHORT_FILE;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
    private static final String FIND_UNKNOWN_SPLICE = "find_unknown_splice";
    private static final String LINE_ELEMENTS_FILE = "line_elements_file";

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
    public final boolean FindUnknownSplice;
    public final String LineElementsFile;

    public FusionCohortConfig(final ConfigBuilder configBuilder)
    {
        MinSampleThreshold = configBuilder.getInteger(FUSION_MIN_SAMPLES);
        MinFragCount = configBuilder.getInteger(FUSION_MIN_FRAGS);
        GenerateCohort = configBuilder.hasFlag(GENERATE_COHORT);
        WriteFilteredFusions = configBuilder.hasValue(WRITE_FILTERED_FUSIONS);
        CompareUnfiltered = configBuilder.hasValue(COMPARE_UNFILTERED);
        RewriteAnnotatedFusions = configBuilder.hasValue(REWRITE_ANNOTATED);
        WriteCombinedFusions = configBuilder.hasValue(WRITE_COMBINED_FUSIONS);
        FindUnknownSplice = configBuilder.hasValue(FIND_UNKNOWN_SPLICE);

        CohortFile = configBuilder.getValue(FUSION_COHORT_FILE);
        ComparisonSource = configBuilder.getValue(FUSION_COMPARISONS);
        LineElementsFile = configBuilder.getValue(LINE_ELEMENTS_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(GENERATE_COHORT, "Generate a fusion cohort file");
        configBuilder.addFlag(WRITE_FILTERED_FUSIONS, "Apply filters to sample fusions and write to file");
        configBuilder.addInteger(FUSION_MIN_SAMPLES, "Min number of samples to support a fusion", 2);
        configBuilder.addInteger(FUSION_MIN_FRAGS, "Min frag count per sample to support a fusion", 2);
        configBuilder.addPath(FUSION_COHORT_FILE, false, "Cohort file previously generated");
        configBuilder.addConfigItem(FUSION_COMPARISONS, "List of sources to compare fusions between");
        configBuilder.addFlag(COMPARE_UNFILTERED, "Included unfiltered fusions in comparison with external tools");
        configBuilder.addFlag(REWRITE_ANNOTATED, "Rewrtite all fusions with cohort and other annotations");
        configBuilder.addFlag(WRITE_COMBINED_FUSIONS, "Write a combined file with passing fusion");
        configBuilder.addFlag(FIND_UNKNOWN_SPLICE, "Find known splice to unknown mappings");
        configBuilder.addPath(LINE_ELEMENTS_FILE, false, "Line elements file");
        addKnownFusionFileOption(configBuilder);
    }
}
