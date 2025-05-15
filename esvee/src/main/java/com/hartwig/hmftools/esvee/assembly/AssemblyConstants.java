package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.common.SvConstants.DEFAULT_MAX_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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
    public static final int MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;
    public static final int DECOY_MAX_MISMATCHES = 3;
    public static final double DECOY_MIN_SCORE_FACTOR = 0.9;
    public static final int ASSEMBLY_MIN_READ_SUPPORT = 2;
    public static final int ASSEMBLY_SPLIT_MIN_READ_SUPPORT = 5;
    public static final double PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC = 0.2;
    public static final int PROXIMATE_REF_SIDE_SOFT_CLIPS = 3;
    public static final int ASSEMBLY_MIN_SOFT_CLIP_LENGTH = MIN_VARIANT_LENGTH;
    public static final int ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH = ASSEMBLY_MIN_SOFT_CLIP_LENGTH / 2;
    public static final int ASSEMBLY_MAX_JUNC_POS_DIFF = 2;
    public static final int ASSEMBLY_REF_READ_MIN_SOFT_CLIP = 10;
    public static final int ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH = 2;
    public static final int ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY = MIN_MAP_QUALITY;
    public static final int ASSEMBLY_MIN_DISTINCT_FRAGS = 2;
    public static final int ASSEMBLY_INDEL_UNLINKED_ASSEMBLY_MIN_LENGTH = 75;

    public static final double DEFAULT_DISC_RATE_INCREMENT = 0.01;
    public static final int DISC_RATE_JUNC_INCREMENT = 1;
    public static final int DISC_RATE_DISC_ONLY_INCREMENT = 3;

    // sequence comparison
    public static final int REPEAT_2_DIFF_COUNT = 11;

    public static final int DEFAULT_ASSEMBLY_MAP_QUAL_THRESHOLD = 10;

    public static final int PRIMARY_ASSEMBLY_MERGE_MISMATCH = 3;
    public static final int PROXIMATE_JUNCTION_DISTANCE = 50;

    // discordant fragment max upper bound is dynamically set from the fragment distribution
    public static int MAX_OBSERVED_CONCORDANT_FRAG_LENGTH = DEFAULT_MAX_CONCORDANT_FRAG_LENGTH;

    // assembly extension
    public static final int ASSEMBLY_READ_OVERLAP_BASES = 20;
    public static final int ASSEMBLY_LINK_OVERLAP_BASES = 50;
    public static final int ASSEMBLY_READ_TRIMMED_OVERLAP_BASES = 30;
    public static final int ASSEMBLY_LINK_DISC_ONLY_OVERLAP_BASES = ASSEMBLY_READ_OVERLAP_BASES;
    public static final int ASSEMBLY_EXTENSION_BASE_MISMATCH = 2;
    public static final int ASSEMBLY_REF_BASE_MAX_GAP = 200;
    public static final int REF_SIDE_MIN_SOFT_CLIP_LENGTH = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

    public static final int LOCAL_ASSEMBLY_MATCH_DISTANCE = 500;
    public static final int MATCH_SUBSEQUENCE_LENGTH = 20;

    // phasing
    public static final int REMOTE_PHASING_MIN_READS = 2;
    public static final int REMOTE_REGION_MERGE_MARGIN = 500;
    public static final double REMOTE_REGION_WEAK_SUPP_PERCENT = 0.1;
    public static final int PHASED_ASSEMBLY_MIN_TI = 30;
    public static final int PHASED_ASSEMBLY_MAX_TI = 1000;
    public static final int PROXIMATE_DEL_LENGTH = 1000;
    public static final int PROXIMATE_DUP_LENGTH = 500;
    public static final double REMOTE_REGION_DISC_READ_BASE_MIN_QUAL_PERC = 0.25;
    public static final int REMOTE_REGION_DISC_READ_BASE_MIN_AS = 75;

    // output
    public static final int DEFAULT_ASSEMBLY_REF_BASE_WRITE_MAX = 200; // for TSV and VCF output, no function impact

    // alignment
    public static final int BWA_PENALTY_ADJUST = 2;
    public static final int ALIGNMENT_MIN_SOFT_CLIP = MIN_VARIANT_LENGTH;
    public static final int ALIGNMENT_MIN_MOD_MAP_QUAL = 10;
    public static final int ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA = 5;
    public static final int ALIGNMENT_CALC_SCORE_FACTOR = 15;
    public static final double ALIGNMENT_CALC_SCORE_THRESHOLD = 0.77;
    public static final int ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH = MIN_ANCHOR_LENGTH;
    public static final int ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH = 50000;
    public static final int ALIGNMENT_LOW_MOD_MQ_QUAL_BOOST = 15;
    public static final int ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH = MIN_ANCHOR_LENGTH;
    public static final int ALIGNMENT_PROXIMATE_DISTANCE = 1000;
    public static final int ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH = 50;
    public static final int ALIGNMENT_RECOVERY_MAX_MD_ERRORS = 4;

    public static final int SHORT_DEL_DUP_INS_LENGTH = 1000;

    // DUX-4 regions
    public static final List<ChrBaseRegion> MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37 = List.of(
            new ChrBaseRegion("4", 190930000, 191030000),
            new ChrBaseRegion("10", 135420000, 135520000 ));

    public static final List<ChrBaseRegion> MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38 = List.of(
            new ChrBaseRegion("chr4", 190060000, 190190000),
            new ChrBaseRegion("chr10", 133660000, 133770000));

    // SSX2 intron 5 region plus SSX2B
    public static final List<ChrBaseRegion> SSX2_REGIONS_V37 = List.of(
            new ChrBaseRegion("X", 52729628, 52731680), // first entry is SSX2
            new ChrBaseRegion("X", 52784877, 52786929));

    public static final List<ChrBaseRegion> SSX2_REGIONS_V38 = List.of(
            new ChrBaseRegion("chrX", 52700578, 52702630),
            new ChrBaseRegion("chrX", 52755800, 52757852));

    public static final Orientation SSX2_GENE_ORIENT = Orientation.REVERSE;
    public static final int SSX2_MAX_MAP_QUAL = 20;

}
