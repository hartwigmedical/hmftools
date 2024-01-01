package com.hartwig.hmftools.esvee;

public final class SvConstants
{
    public static final String APP_NAME = "Esvee";

    // file related
    public static final String ASSEMBLY_BAM_FILE_ID = ".assembly.bam";
    public static final String SV_PREP_JUNCTIONS_FILE_ID = ".sv_prep.junctions.tsv";
    public static final String DEFAULT_HTML_SUMMARY_DIR = "html";
    public static final int MAX_HTML_SUMMARIES = 10000;

    public static final String BAM_HEADER_SAMPLE_ID_TAG = "sampleId";

    public static final int BAM_READ_JUNCTION_BUFFER = 1000;

    public static final int MAX_DUP_LENGTH = 6; // then classified as an INS

    // When a base is considered in the context of a single read, at what quality level do we start to see this base as low-quality"
    public static int LOW_BASE_QUAL_THRESHOLD = 26; // lowQualBaseThreshold

    // When a base has more than 1 piece of supporting evidence, what is the cumulative total below which the base is low quality anyway")
    public static int LOW_BASE_QUAL_CUMULATIVE_THRESHOLD = 39; // was lowBaseQualCumulativeThreshold

    // The average base qual below which we consider the entire read (or section of a read) to be low quality
    public static int AVG_BASE_QUAL_THRESHOLD = 30; // averageQualityThreshold

    public static final int MAX_MISMATCH_STRONG_SUPPORT = 1;

    public static final int MAX_MISMATCH_WEAK_SUPPORT = 4;

    public static final int MAX_MISMATCH_DEDUPING_ASSEMBLIES = 5;

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

    public static final int MIN_MAPQ_START_JUNCTION = 20;

    // When trimming polgG/C, how many Gs/Cs must appear on the edge of a read to be trimmed
    public static final int NORMALISERPOLYGLENGTH = 4;

    // Indels this size or larger near the edge of a read will be converted to soft-clips
    public static final int NORMALISERINDELMINSIZETOSOFTCLIP = 6;

    // How close to the edge is \"near\" for an indels
    public static final int NORMALISERINDELMAXEDGEDISTANCE = 16;

    // Not a MapQ. Ignore any alignments during extension that come from BWA with a score less than this
    public static final int ALIGNERMINSCORE = 20;

    // If there is more than one candidate for extending an assembly alignment, ignore any that insert this many more bases than the best
    public static final int ALIGNEREXTENSIONINSERTTOLERANCE = 16;

    // We don't attempt to call the aligner if we have less than this many bases
    public static final int ALIGNERMINBASES = 20;

    // When extending alignments, if we have a candidate match within this many bases of the existing neighbour, prioritise that alignment
    public static final int ALIGNERMAXDISTANCETOCONSIDERNEARBY = 2000;

    // What is the smallest insertion/deletion size to call
    public static final int CALLERMINSIZETOCALL = 32;

    // The threshold below which a LOW_OVERHANG filter will be applied to the VCF
    public static final int VCFLOWOVERHANGTHRESHOLD = 30;

    // The threshold below which a LOW_QUALITY filter will be applied to the VCF
    public static final int VCFLOWQUALITYTHRESHOLD = 40;


    public static final boolean EXTEND_PRIMARIES = false;

    // performance and logging related
    public static final int TASK_LOG_COUNT = 100;

    // Whether individual operations are timed to prevent slow processing
    public static final boolean TIMEOUTS_ENABLED = false;
    public static final int PRIMARY_TIMEOUT = 5000;
    public static final int EXTENSION_TIMEOUT = 5000;


}
