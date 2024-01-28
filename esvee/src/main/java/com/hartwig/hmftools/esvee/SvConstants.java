package com.hartwig.hmftools.esvee;

public final class SvConstants
{
    public static final String APP_NAME = "Esvee";

    // file related
    public static final String SV_PREP_JUNCTIONS_FILE_ID = ".sv_prep.junctions.tsv";
    public static final String DEFAULT_HTML_SUMMARY_DIR = "html";
    public static final String REF_GENOME_IMAGE_EXTENSION = ".img";
    public static final int MAX_HTML_SUMMARIES = 10000;

    public static final String BAM_HEADER_SAMPLE_ID_TAG = "sampleId";

    // BAM reading
    public static final int BAM_READ_JUNCTION_BUFFER = 1000;

    // read adjustments
    public static final int POLY_G_TRIM_LENGTH = 4;
    public static final int INDEL_TO_SC_MIN_SIZE_SOFTCLIP = 6;
    public static final int INDEL_TO_SC_MAX_EDGE_DISTANCE = 16;

    // primary assembly
    public static final int READ_SOFT_CLIP_JUNCTION_BUFFER = 2;

    public static final int READ_FILTER_MIN_JUNCTION_MAPQ = 20;
    public static final int READ_FILTER_MIN_ALIGNED_BASES = 30; // new (ie previously hard-coded)

    public static final int PRIMARY_ASSEMBLY_MIN_LENGTH = 10;
    public static final int PRIMARY_ASSEMBLY_MIN_READ_SUPPORT = 2; // new
    public static final int PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL = 60; // new

    // primary assembly deduplication
    public static final int PRIMARY_ASSEMBLY_READ_MAX_BASE_MISMATCH = 1;
    public static final int PRIMARY_ASSEMBLY_MERGE_MISMATCH = 3;
    public static final int PRIMARY_ASSEMBLY_MERGE_READ_OVERLAP = 2;

    public static final int PROXIMATE_JUNCTION_DISTANCE = 10;
    public static final int PROXIMATE_JUNCTION_OVERLAP = 100;

    // filters
    public static final int PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH = 32;

    // common
    public static int AVG_BASE_QUAL_THRESHOLD = 30; // original name: averageQualityThreshold
    public static int LOW_BASE_QUAL_THRESHOLD = 26; // original name: lowQualBaseThreshold



    public static final int SUPPORT_MAX_MISMATCH_STRONG = 1;
    public static final int SUPPORT_MAX_MISMATCH_WEAK = 4;
    public static final int SUPPORT_MAX_MISMATCH_DEDUPING_ASSEMBLIES = 5;


    // variant calling
    public static final int MAX_DUP_LENGTH = 6; // then classified as an INS


    // OLD CODE - may still be used
    // read support
    public static final int PRIMARY_ASSEMBLY_WEAK_SUPPORT_MIN_BASES = 3; // new


    // UNCLASSIFIED


    // When a base is considered in the context of a single read, at what quality level do we start to see this base as low-quality"

    // When a base has more than 1 piece of supporting evidence, what is the cumulative total below which the base is low quality anyway")
    public static int LOW_BASE_QUAL_CUMULATIVE_THRESHOLD = 39; // was lowBaseQualCumulativeThreshold

    // The average base qual below which we consider the entire read (or section of a read) to be low quality



    public static final int MAX_DISTANCE_EDUPE_ASSEMBLIES = 50;

    public static final boolean TRY_EXTENDING_USING_DISCORDANT_READS = true;

    public static final int DISCORDANT_FRAGMENT_LENGTH = 1000;

    public static final int DISCORDANT_PAIR_SEARCH_DISTANCE = 400;

    public static final int DISCORDANT_PAIR_MIN_MAPQ = 31;

    public static final int ASSEMBLY_EXTENSION_MIN_MATCH_BASES = 20;

    public static final int ASSEMBLY_EXTENSION_MAX_REPEAT_SCORE = 4;

    public static final int ASSEMBLY_EXTENSION_MAX_MISMATCH = 1;


    public static final int ASSEMBLY_EXTENSION_MIN_SUPPORT = 8;
    public static final int ASSEMBLY_EXTENSION_MIN_SUPPORT_FINAL = 4;



    public static final int MAX_MISMATCH_FOLDING = 1;

    // If a variant or assembly has less than this many fragments support it, it is dropped
    public static final int MIN_READS_SUPPORT_ASSEMBLY = 2;



    // Not a MapQ. Ignore any alignments during extension that come from BWA with a score less than this
    public static final int ALIGNER_MIN_SCORE = 20;

    // If there is more than one candidate for extending an assembly alignment, ignore any that insert this many more bases than the best
    public static final int ALIGNER_EXTENSION_INSERT_TOLERANCE = 16;

    // We don't attempt to call the aligner if we have less than this many bases
    public static final int ALIGNER_MIN_BASES = 20;

    // When extending alignments, if we have a candidate match within this many bases of the existing neighbour, prioritise that alignment
    public static final int ALIGNER_NEARBY_MAX_DISTANCE = 2000;

    // What is the smallest insertion/deletion size to call
    public static final int VARIANT_MIN_LENGTH = 32;

    // The threshold below which a LOW_OVERHANG filter will be applied to the VCF
    public static final int LOW_OVERHANG_THRESHOLD = 30;

    // The threshold below which a LOW_QUALITY filter will be applied to the VCF
    // FIXME: make config
    public static final int LOW_QUALITY_THRESHOLD = 40;


    public static final boolean EXTEND_PRIMARIES = false;

    // performance and logging related
    public static final int TASK_LOG_COUNT = 100;

    // Whether individual operations are timed to prevent slow processing
    public static final boolean TIMEOUTS_ENABLED = false;
    public static final int PRIMARY_TIMEOUT = 5000;
    public static final int EXTENSION_TIMEOUT = 5000;


}
