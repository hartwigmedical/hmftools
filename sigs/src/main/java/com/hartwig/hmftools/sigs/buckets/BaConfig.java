package com.hartwig.hmftools.sigs.buckets;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class BaConfig
{
    // input files
    public static final String BA_EXT_SAMPLE_DATA_FILE = "ba_ext_data_file";
    public static final String BA_SAMPLE_CALC_DATA_FILE = "ba_sam_calc_data_file";
    public static final String BA_PREDEFINED_SIGS = "ba_predefined_sigs_file";

    public static final double MAX_ELEVATED_PROB = 1e-12;

    public static final double MIN_DISCOVERY_SAMPLE_COUNT = 0.0001; // % of cohort to consider a pair of samples similar (5K at 0.01% of 55M)
    public static final int MIN_CANCER_TYPE_SAMPLES = 10;

    public static double MIN_GROUP_ALLOC_PERCENT = 0.10; // required elevated count % allocation of sample during discovery
    public static final String BA_MIN_GROUP_ALLOC_PERCENT = "ba_min_grp_alloc_pct";
    public static double MIN_GROUP_ALLOC_PERCENT_LOWER = 0.05; // required total count % allocation of sample during final fit (background sig excluded)
    public static final String BA_MIN_GROUP_ALLOC_PERCENT_LOWER = "ba_min_grp_alloc_pct_low";

    public static final double SKIP_ALLOC_FACTOR = 2.0; // skip adding a sample to a group if another candidate not chosen to allocate X times as much

    // configuration state
    public final double HighCssThreshold; // CSS level for samples or groups to be consider similar
    private static final  String BA_CSS_HIGH_THRESHOLD = "ba_css_high";
    public static final double DEFAULT_CSS_THRESHOLD = 0.995;

    public final int MutationalLoadCap; // defines the maximum for background counts
    private static final  String BA_MUT_LOAD_CAP = "ba_mut_load_cap";

    public final int MinBucketCountOverlap; // used for pairings of samples with reduced bucket overlap
    private static final  String BA_MIN_BUCKET_COUNT_OVERLAP = "ba_min_bc_overlap";

    public final int MinSampleAllocCount; // hard lower limit to allocate a sample to a group
    private static final String BA_MIN_SAM_ALLOC_COUNT = "ba_min_sample_alloc_count";

    public final boolean ApplyNoise; // whether to factor Poisson noise into the sample counts and fits, default = true
    public static final double PERMITTED_PROB_NOISE = 1e-5;
    public static final double MAX_NOISE_TO_SAMPLE_RATIO = 5; // per sample, the sum of potential (unallocated) bucket noise counts cannot exceeds this multiple of variant total
    public static double MAX_NOISE_ALLOC_PERCENT = 0.20; // per sample, the total which can be allocated to noise vs total variants
    public static final String BA_MAX_NOISE_ALLOC_PERCENT = "ba_max_noise_pct";

    public final boolean UseRatioRanges; // default false
    private static final  String BA_USE_RATIO_RANGES = "ba_use_ratio_ranges";
    public final double RatioRangeLimit;
    private static final  String BA_RATIO_RANGE = "ba_ratio_range";
    public static final double DEFAULT_SIG_RATIO_RANGE_PERCENT = 0.1; // if ratio ranges are used, this percent width can be applied

    public final double MutLoadWeightFactor; // how much weight to give a sample's counts towards bucket ratios (1 = full)
    private static final  String BA_MUT_LOAD_WEIGHT_FACTOR = "ba_mut_load_wf";

    public static final double UNIQUE_SIG_CSS_THRESHOLD = 0.8; // for unique sig discovery, sample's bucket ratios must be no more similar than this to any other sig
    public static final double UNIQUE_SIG_MIN_ALLOC_PERCENT = 0.3; // for unique sig discovery, the unique sig must allocation at least this % to the new sig

    public static final double MAJOR_GROUP_ALLOC_PERC = 0.02; // define a 'MAJOR' sig as contributing at least this level to the whole cohort
    public static final double MAJOR_GROUP_SAMPLE_PERC = 0.05; // in total count and percentage terms

    public static final double SIG_SIMILAR_CSS = 0.90; // for sig and bucket group similarity
    public static final double DOMINANT_CATEGORY_PERCENT = 0.7; // mark a group with a category if X% of samples in it have this attribute (eg cancer type, UV)

    public static final double SAMPLE_ALLOCATED_PERCENT = 0.995; // maximum to attempt to allocate a sample to sigs

    public final boolean UseBackgroundCounts; // whether to make a distinction between background and elevated counts
    private static final String BA_USE_BACKGROUND_SIGS = "ba_use_background_sigs"; // true by default
    public final boolean RunBackgroundAnalysis;
    private static final  String BA_RUN_BACKGROUND_ANALYSIS = "ba_run_background_analysis"; // false by default

    public final int RunCount; //  number of iterations searching for potential bucket groups / sigs
    private static final  String BA_RUN_COUNT = "ba_run_count";

    public final int ExcessDiscoveryRunId; // run iteration from which to begin discovery using excess-unalloc method
    private static final  String BA_EXCESS_GRP_RUN_INDEX = "ba_excess_grp_run_index";

    public final boolean ThrottleDiscovery; // whether to allow discovery after each run
    public final int UniqueDiscoveryRunId;

    public final double GroupMergeScoreThreshold; // after discovery is complete, attempt to merge similar sigs if score exceeds this
    private static final  String BA_MERGE_SIG_SCORE = "ba_merge_sig_score";

    public final  int MaxProposedSigs;
    private static final  String BA_MAX_PROPOSED_SIGS = "ba_max_proposed_sigs";

    public final  int ApplyPredefinedSigCount; // how many loaded sigs to apply prior to discovery (optimisation)
    public static final String BA_PREDEFINED_SIG_APPLY_COUNT = "ba_predefined_sig_apply_count";

    public static final String CANCER_TYPE_OTHER = "Other";
    // data filters
    public final String SpecificCancer;
    private static final  String BA_SPECIFIC_CANCER = "ba_specific_cancer";

    public final String MsiFilter;
    private static final String BA_MSI_FILTER = "ba_msi_filter";

    // performance related
    public static final int MAX_CANDIDATE_GROUPS = 1500; // in place for speed and memory considerations

    // logging and debug
    public final List<Integer> SampleWatchList;
    private static final  String BA_LOG_SPEC_SAMPLES = "ba_log_sample_ids";

    public boolean LogVerbose;

    public BaConfig(final CommandLine cmd)
    {
        LogVerbose = false;

        UseBackgroundCounts = Boolean.parseBoolean(cmd.getOptionValue(BA_USE_BACKGROUND_SIGS, "true"));
        RunBackgroundAnalysis = cmd.hasOption(BA_RUN_BACKGROUND_ANALYSIS);

        RunCount = Integer.parseInt(cmd.getOptionValue(BA_RUN_COUNT, "25"));
        ExcessDiscoveryRunId = Integer.parseInt(cmd.getOptionValue(BA_EXCESS_GRP_RUN_INDEX, "-1"));
        UniqueDiscoveryRunId = -1;
        ThrottleDiscovery = false;

        MaxProposedSigs = Integer.parseInt(cmd.getOptionValue(BA_MAX_PROPOSED_SIGS, "0"));
        ApplyPredefinedSigCount = Integer.parseInt(cmd.getOptionValue(BA_PREDEFINED_SIG_APPLY_COUNT, "0"));
        UseRatioRanges = cmd.hasOption(BA_USE_RATIO_RANGES);

        MutLoadWeightFactor = Double.parseDouble(cmd.getOptionValue(BA_MUT_LOAD_WEIGHT_FACTOR, "1"));

        MinSampleAllocCount = Integer.parseInt(cmd.getOptionValue(BA_MIN_SAM_ALLOC_COUNT, "1"));
        MinBucketCountOverlap = Integer.parseInt(cmd.getOptionValue(BA_MIN_BUCKET_COUNT_OVERLAP, "3"));
        HighCssThreshold = Double.parseDouble(cmd.getOptionValue(BA_CSS_HIGH_THRESHOLD, String.valueOf(DEFAULT_CSS_THRESHOLD)));
        MutationalLoadCap = Integer.parseInt(cmd.getOptionValue(BA_MUT_LOAD_CAP, "10000"));
        RatioRangeLimit = Double.parseDouble(cmd.getOptionValue(BA_RATIO_RANGE, String.valueOf(DEFAULT_SIG_RATIO_RANGE_PERCENT)));

        MAX_NOISE_ALLOC_PERCENT = Double.parseDouble(cmd.getOptionValue(BA_MAX_NOISE_ALLOC_PERCENT, String.valueOf(MAX_NOISE_ALLOC_PERCENT)));

        ApplyNoise = MAX_NOISE_ALLOC_PERCENT > 0;

        MIN_GROUP_ALLOC_PERCENT = Double.parseDouble(cmd.getOptionValue(BA_MIN_GROUP_ALLOC_PERCENT, String.valueOf(MIN_GROUP_ALLOC_PERCENT)));
        MIN_GROUP_ALLOC_PERCENT_LOWER = Double.parseDouble(cmd.getOptionValue(BA_MIN_GROUP_ALLOC_PERCENT_LOWER, String.valueOf(MIN_GROUP_ALLOC_PERCENT_LOWER)));

        GroupMergeScoreThreshold = Double.parseDouble(cmd.getOptionValue(BA_MERGE_SIG_SCORE, "0"));

        SpecificCancer = cmd.getOptionValue(BA_SPECIFIC_CANCER, "");
        MsiFilter = cmd.getOptionValue(BA_MSI_FILTER, "");

        SampleWatchList = Lists.newArrayList();

        if(cmd.hasOption(BA_LOG_SPEC_SAMPLES))
        {
            SampleWatchList.addAll(Arrays.stream(cmd.getOptionValue(BA_LOG_SPEC_SAMPLES).split(";"))
                    .map(x -> Integer.parseInt(x))
                    .collect(Collectors.toList()));
        }
    }

    public boolean logSample(int sampleId)
    {
        if(SampleWatchList.contains(sampleId))
            return true;

        return false;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(BA_USE_BACKGROUND_SIGS, true, "Default true. Whether to calculate and use background sigs");
        options.addOption(BA_RUN_BACKGROUND_ANALYSIS, false, "Analysis background sigs by cancer type");
        options.addOption(BA_EXT_SAMPLE_DATA_FILE, true, "Sample external data");
        options.addOption(BA_CSS_HIGH_THRESHOLD, true, "Cosine sim for high-match test");
        options.addOption(BA_MAX_PROPOSED_SIGS, true, "Maximum number of bucket groups to turn into proposed sigs");
        options.addOption(BA_SAMPLE_CALC_DATA_FILE, true, "Optional: file containing computed data per sample");
        options.addOption(BA_SPECIFIC_CANCER, true, "Optional: Only process this cancer type");
        options.addOption(BA_RUN_COUNT, true, "Number of search iterations");
        options.addOption(BA_PREDEFINED_SIGS, true, "Predefined sigs to use during discovery");
        options.addOption(BA_PREDEFINED_SIG_APPLY_COUNT, true, "How many predefined sigs to apply prior to discovery");
        options.addOption(BA_MIN_SAM_ALLOC_COUNT, true, "Min count to allocate a sample to a group");
        options.addOption(BA_MUT_LOAD_CAP, true, "Mutational load cap used in background counts calc");
        options.addOption(BA_EXCESS_GRP_RUN_INDEX, true, "Run id for excess-unalloc group logic to kick in");
        options.addOption(BA_MIN_BUCKET_COUNT_OVERLAP, true, "Min buckets for candidate group discovery");
        options.addOption(BA_USE_RATIO_RANGES, false, "Allow a computed range around sig ratios");
        options.addOption(BA_RATIO_RANGE, true, "Ratio range limit");
        options.addOption(BA_MAX_NOISE_ALLOC_PERCENT, true, "Max noise per bucket as % of count");
        options.addOption(BA_MUT_LOAD_WEIGHT_FACTOR, true, "Dampen higher mutational load samples (1 = no dampening)");
        options.addOption(BA_MSI_FILTER, true, "Use 'Include' to only look at MSI samples, or 'Exclude' to exclude them");
        options.addOption(BA_MERGE_SIG_SCORE, true, "After discovery, merge similar sigs before final fit");
        options.addOption(BA_MIN_GROUP_ALLOC_PERCENT, true, "Only allocate to a group if exceeds this percent (background excepted)");
        options.addOption(BA_MIN_GROUP_ALLOC_PERCENT_LOWER, true, "Lower threshold for group allocation");
        options.addOption(BA_LOG_SPEC_SAMPLES, true, "Set of samples IDs to log, separated by ';'");
    }

}
