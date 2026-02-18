package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion

object CiderConstants
{
    // for partially rearranged VDJs, we designate that the "interesting"
    // sequence is 60 bases before the J anchor or 60 bases after the V anchor
    const val PARTIAL_VDJ_UNANCHORED_LENGTH_BASES: Int = 60

    // following are constants for BLOSUM search
    const val MAX_BLOSUM_DIFF_PER_AA: Int = 3
    const val BLOSUM_SIMILARITY_SCORE_CONSTANT: Int = 6
    const val BLOSUM_UNKNOWN_BASE_PENALTY: Int = 2

    const val CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 10
    const val VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 3

    // The minimum number of bases a read must overlap the CDR3 sequence (i.e. sequence past the anchor) to be used.
    const val MIN_READ_CDR3_OVERLAP: Int = 20

    // minimum amount of bases a read must overlap with a layout for it to be added
    const val LAYOUT_MIN_READ_OVERLAP_BASES: Int = 20
    const val MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES: Int = 20

    // Maximum number of mismatches between the consensus and a mate read attempting to extend it.
    const val LAYOUT_MATE_EXTEND_MISMATCHES_MAX = 1

    // word size of 8 is experimentally decided to provide good speed up in finding VJ overlaps
    const val VJ_JOIN_HASH_WORD_SIZE: Int = 8
    const val MIN_ANCHOR_LENGTH_BASES: Int = 3
    const val LAYOUT_MIN_SUPPORT_TO_SEAL_NODE: Int = 5

    // minimum fraction of reads that support a certain base in a layout
    // this is very important such that we do not let one or two reads stop
    // layouts from being joined together
    const val MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION: Double = 0.02
    const val MIN_POLY_G_TRIM_COUNT: Int = 2
    const val POLY_G_TRIM_EXTRA_BASE_COUNT: Int = 5

    const val MIN_NON_SPLIT_READ_STRADDLE_LENGTH: Int = 30

    // From some testing, 50k query sequences used 10GB of memory.
    // Want to limit it to about use only a few GB.
    const val BWAMEM_BATCH_SIZE = 20000

    // Require a match of minimum ~20 bases. If we want to match D segment that is shorter
    // we will need a higher cut off, maybe 10, but will get many false positive hits that are longer but more mismatches
    const val ANNOTATION_ALIGN_SCORE_MIN = 19
    const val ANNOTATION_VDJ_FLANK_BASES = 300
    // filter out matches that have too low identity
    // reason for doing this is that we use match/mismatch of 1/-4, in worst case we can
    // get 1 mismatch for every 4 matches, and could find alignments with 80% identity.
    // those are probably too different to use. We use 90% for V / J identities, and 95%
    // cut off for full match
    const val ANNOTATION_VJ_IDENTITY_MIN = 0.9
    const val ANNOTATION_MATCH_REF_IDENTITY = 0.95
    // We require that the alignment covers the anchor boundary base, +/- this many bases.
    // This is because sometimes a valid rearrangement can chop off the last few bases of the anchor.
    const val ANNOTATION_ANCHOR_BASE_TOLERANCE = 10

    // For the scoring function, the match/mismatch score 1/-4 optimizes the scoring for 100% identical sequences and 1/-1 for 75% identical
    // sequences. The default for NCBI Blastn is 2/-3, which is optimal for 89% identical sequences. BWA uses 1/-4.
    // There is also gap opening and gap extension. BWA uses gap opening of -6 and gap extension of -1.
    // For blastn default scoring, see: https://www.ncbi.nlm.nih.gov/books/NBK279684/
    // we set it to 1/-4/-5/-2 which is optimal for 100% identical sequences
    // My test shows that this scoring would mostly prefer shorter matches with higher identity than longer matches with lower identity.
    const val ALIGNER_MATCH_SCORE = 1
    const val ALIGNER_MISMATCH_SCORE = -4
    const val ALIGNER_GAP_OPENING_SCORE = -5
    const val ALIGNER_GAP_EXTEND_SCORE = -2

    const val ALIGNER_WORD_SIZE = 9

    // blast uses v38
    val BLAST_REF_GENOME_VERSION = RefGenomeVersion.V38
    // This limits the number of hit each query can get. Necessary to protect against edge cases
    const val BLASTN_MAX_TARGET_SEQUENCES = 5000
    val BLASTN_CHROMOSOME_ASSEMBLY_REGEX = Regex("^Homo sapiens chromosome (\\w+).*, GRCh38.p13 (.+)$")
    val BLASTN_PRIMARY_ASSEMBLY_NAME = "Primary Assembly".intern()

    // Amino acids sequences which are known to match the reference genome but are not detected by alignment.
    // This exists because of switching from Blastn to BWA-MEM, and some sequences cause significant discrepancies.
    // These sequences are not in the ref genome for BWA-MEM to align to.
    val MATCHES_REF_KNOWN_CDR3_AA = listOf(
        Regex("CTXGPKXELRT.*")
    )
}