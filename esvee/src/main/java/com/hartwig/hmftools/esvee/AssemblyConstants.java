package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.esvee.common.SvConstants.DEFAULT_DISCORDANT_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

public final class AssemblyConstants
{
    // BAM reading
    public static final int BAM_READ_JUNCTION_BUFFER = 1000;

    // read adjustments
    public static final int INDEL_TO_SC_MIN_SIZE_SOFTCLIP = MIN_INDEL_SUPPORT_LENGTH;
    public static final int INDEL_TO_SC_MAX_SIZE_SOFTCLIP = MIN_VARIANT_LENGTH - 1;
    public static final int POLY_G_TRIM_LENGTH = 4;
    public static final double LOW_BASE_TRIM_PERC = 0.35;
    public static final int UNMAPPED_TRIM_THRESHOLD = 40;

    // primary assembly
    public static final int READ_SOFT_CLIP_JUNCTION_BUFFER = 2;
    public static final int MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;;
    public static final int DECOY_MAX_MISMATCHES = 3;
    public static final int ASSEMBLY_MIN_READ_SUPPORT = 2;
    public static final int ASSEMBLY_SPLIT_MIN_READ_SUPPORT = 5;
    public static final double PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC = 0.2;
    public static final int PROXIMATE_REF_SIDE_SOFT_CLIPS = 3;
    public static final int ASSEMBLY_MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;
    public static final int ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH = ASSEMBLY_MIN_SOFT_CLIP_LENGTH / 2;
    public static final int ASSEMBLY_MAX_JUNC_POS_DIFF = 2;
    public static final int ASSEMBLY_REF_READ_MIN_SOFT_CLIP = 10;
    public static final int ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH = 2;

    public static final int DEFAULT_ASSEMBLY_MAP_QUAL_THRESHOLD = 10;

    public static final int PRIMARY_ASSEMBLY_MERGE_MISMATCH = 3;
    public static final int PROXIMATE_JUNCTION_DISTANCE = 50;

    // discordant fragment max upper bound is dynamically set from the fragment distribution
    public static int DISCORDANT_FRAGMENT_LENGTH = DEFAULT_DISCORDANT_FRAGMENT_LENGTH;

    // assembly extension
    public static final int ASSEMBLY_READ_OVERLAP_BASES = 20;
    public static final int ASSEMBLY_LINK_OVERLAP_BASES = 50;
    public static final int ASSEMBLY_EXTENSION_BASE_MISMATCH = 2;
    public static final int ASSEMBLY_REF_BASE_MAX_GAP = 200;
    public static final int REF_SIDE_MIN_SOFT_CLIP_LENGTH = MIN_SOFT_CLIP_LENGTH;

    public static final int LOCAL_ASSEMBLY_MATCH_DISTANCE = 500;
    public static final int MATCH_SUBSEQUENCE_LENGTH = 20;

    // phasing
    public static final int REMOTE_PHASING_MIN_READS = 2;
    public static final int REMOTE_REGION_MERGE_MARGIN = 100;
    public static final int REMOTE_REGION_REF_MIN_READS = REMOTE_PHASING_MIN_READS;
    public static final double REMOTE_REGION_REF_MIN_READ_PERCENT = 0.1;
    public static final double REMOTE_REGION_WEAK_SUPP_PERCENT = 0.1;
    public static final int PHASED_ASSEMBLY_MIN_TI = 30;
    public static final int PHASED_ASSEMBLY_MAX_TI = 1000;
    public static final int PROXIMATE_DEL_LENGTH = 1000;
    public static final int PROXIMATE_DUP_LENGTH = 500;

    // output
    public static final int DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX = 200; // for TSV and VCF output, no function impact

    // alignment
    public static final int ALIGNMENT_MIN_SOFT_CLIP = 30;
    public static final int ALIGNMENT_MIN_MOD_MAP_QUAL = 10;
    public static final int ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA = 5;
    public static final int ALIGNMENT_CALC_SCORE_FACTOR = 15;
    public static final double ALIGNMENT_CALC_SCORE_THRESHOLD = 0.85;
    public static final int ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH = MIN_ANCHOR_LENGTH;
    public static final int ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH = 50000;
    public static final int ALIGNMENT_LOW_MOD_MQ_QUAL_BOOST = 15;
    public static final int ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH = MIN_ANCHOR_LENGTH;

    public static final int SHORT_DEL_DUP_INS_LENGTH = 1000;
}
