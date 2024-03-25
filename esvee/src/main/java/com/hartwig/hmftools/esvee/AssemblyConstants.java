package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

public final class AssemblyConstants
{
    public static final String APP_NAME = "Esvee";

    // file related
    public static final String REF_GENOME_IMAGE_EXTENSION = ".img";

    public static final String BAM_HEADER_SAMPLE_ID_TAG = "sampleId";

    // BAM reading
    public static final int BAM_READ_JUNCTION_BUFFER = 1000;

    // read adjustments
    public static final int INDEL_TO_SC_MIN_SIZE_SOFTCLIP = 10;
    public static final int INDEL_TO_SC_MAX_SIZE_SOFTCLIP = MIN_VARIANT_LENGTH - 1;

    // primary assembly
    public static final int READ_SOFT_CLIP_JUNCTION_BUFFER = 2;
    public static final int MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;;

    public static final int PROXIMATE_DEL_LENGTH = 1000;
    public static final int PROXIMATE_DUP_LENGTH = 500;

    public static final int DECOY_MAX_MISMATCHES = 3;

    public static final int PRIMARY_ASSEMBLY_MIN_READ_SUPPORT = 2;
    public static final int PROXIMATE_REF_SIDE_SOFT_CLIPS = 3;
    public static final int PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;

    // primary assembly deduplication
    public static final int PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH = 1;
    public static final int PRIMARY_ASSEMBLY_SUPPORT_MISMATCH = 2;
    public static final int PRIMARY_ASSEMBLY_MERGE_MISMATCH = 3;

    public static final int PROXIMATE_JUNCTION_DISTANCE = 50;

    public static final double REMOTE_REGION_WEAK_SUPP_PERCENT = 0.1;

    // assembly extension
    public static final int ASSEMBLY_EXTENSION_OVERLAP_BASES = 20;
    public static final int ASSEMBLY_EXTENSION_BASE_MISMATCH = 2;
    public static final int REF_SIDE_MIN_SOFT_CLIP_LENGTH = MIN_SOFT_CLIP_LENGTH;

    // phasing
    public static final int REMOTE_PHASING_MIN_READS = 2;
    public static final int PHASED_ASSEMBLY_OVERLAP_BASES = 100;
    public static final int PHASED_ASSEMBLY_JUNCTION_OVERLAP = 50;
    public static final int PHASED_ASSEMBLY_MAX_TI = 1000;

    public static final int DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX = 200; // for TSV and VCF output, no function impact

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
