package com.hartwig.hmftools.esvee;

public final class SvConstants
{
    public static final String APP_NAME = "Esvee";

    // file related
    public static final String SV_PREP_JUNCTIONS_FILE_ID = ".sv_prep.junctions.tsv";
    public static final String REF_GENOME_IMAGE_EXTENSION = ".img";

    public static final String BAM_HEADER_SAMPLE_ID_TAG = "sampleId";

    // BAM reading
    public static final int BAM_READ_JUNCTION_BUFFER = 1000;

    // common
    public static final int MIN_VARIANT_LENGTH = 32;

    // read adjustments
    public static final int POLY_G_TRIM_LENGTH = 4;
    public static final int INDEL_TO_SC_MIN_SIZE_SOFTCLIP = 6;
    public static final int INDEL_TO_SC_MAX_SIZE_SOFTCLIP = MIN_VARIANT_LENGTH - 1;
    public static final int INDEL_TO_SC_MAX_EDGE_DISTANCE = 50;
    public static final double LOW_BASE_TRIM_PERC = 0.3;

    // primary assembly
    public static final int READ_SOFT_CLIP_JUNCTION_BUFFER = 2;

    public static final int PROXIMATE_DEL_LENGTH = 1000;
    public static final int PROXIMATE_DUP_LENGTH = 500;

    public static final int DECOY_MAX_MISMATCHES = 3;

    public static final int PRIMARY_ASSEMBLY_MIN_LENGTH = 10;
    public static final int PRIMARY_ASSEMBLY_MIN_READ_SUPPORT = 2; // new
    public static final int PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL = 60; // new
    public static final int PROXIMATE_REF_SIDE_SOFT_CLIPS = 3;
    public static final int PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;

    // primary assembly deduplication
    public static final int PRIMARY_ASSEMBLY_MAX_BASE_MISMATCH = 1;
    public static final int PRIMARY_ASSEMBLY_MERGE_MISMATCH = 3;
    public static final int PRIMARY_ASSEMBLY_MERGE_READ_OVERLAP = 2;

    public static final int PROXIMATE_JUNCTION_DISTANCE = 50;
    public static final int PROXIMATE_JUNCTION_OVERLAP = 100;

    // filters
    public static int AVG_BASE_QUAL_THRESHOLD = 30;
    public static int LOW_BASE_QUAL_THRESHOLD = 26;

    public static final double REMOTE_REGION_WEAK_SUPP_PERCENT = 0.1;

    // assembly extension
    public static final int ASSEMBLY_EXTENSION_OVERLAP_BASES = 20;
    public static final int ASSEMBLY_EXTENSION_BASE_MISMATCH = 2;

    // phased assembly overlaps
    public static final int PHASED_ASSEMBLY_OVERLAP_BASES = 100;
    public static final int PHASED_ASSEMBLY_JUNCTION_OVERLAP = 50;
    public static final int PHASED_ASSEMBLY_MAX_TI = 1000;

    // CHECK: share with SvPrep, possibly shorten
    public static final int MIN_INDEL_SUPPORT_LENGTH = 15;
    public static final int MIN_INDEL_LENGTH = MIN_VARIANT_LENGTH;

    // variant calling
    public static final int MAX_DUP_LENGTH = 6; // then classified as an INS

    public static final int DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX = 200; // for TSV and VCF output, no function impact


    public static final int DISCORDANT_FRAGMENT_LENGTH = 1000;

    // alignment - to check

    // Not a MapQ. Ignore any alignments during extension that come from BWA with a score less than this
    public static final int ALIGNER_MIN_SCORE = 20;

    // If there is more than one candidate for extending an assembly alignment, ignore any that insert this many more bases than the best
    public static final int ALIGNER_EXTENSION_INSERT_TOLERANCE = 16;

    // We don't attempt to call the aligner if we have less than this many bases
    public static final int ALIGNER_MIN_BASES = 20;

    // When extending alignments, if we have a candidate match within this many bases of the existing neighbour, prioritise that alignment
    public static final int ALIGNER_NEARBY_MAX_DISTANCE = 2000;

    // The threshold below which a LOW_OVERHANG filter will be applied to the VCF
    public static final int LOW_OVERHANG_THRESHOLD = 30;

    // The threshold below which a LOW_QUALITY filter will be applied to the VCF
    public static final int LOW_QUALITY_THRESHOLD = 40;
}
