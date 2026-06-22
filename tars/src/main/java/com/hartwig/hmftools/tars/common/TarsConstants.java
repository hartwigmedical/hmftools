package com.hartwig.hmftools.tars.common;

// All Tars constants in one place: tool identity, output file id, the transcript-contig naming scheme, and the
// algorithm/policy thresholds every liftback pass shares. Config-key strings stay with the config class; local
// plumbing constants stay with their owning class.
public final class TarsConstants
{
    // Tool identity.
    public static final String APP_NAME = "Tars";

    // Output file id, e.g. <sample>.tars.summary.tsv (mirrors Redux's FILE_ID).
    public static final String FILE_ID = "tars";

    // Transcript-contig naming. An alt contig is the chromosome name plus the _tx suffix (e.g. chr1_tx); the
    // delimiter and prefix describe the per-transcript naming used when building the contig FASTA.
    public static final String CONTIG_NAME_DELIM = "_";
    public static final String CONTIG_NAME_PREFIX = "ens";
    public static final String ALT_CONTIG_SUFFIX = "_tx";

    // SpliceFastaBuilder output file ids.
    public static final String TRANSCRIPT_CONTIGS_FILE_ID = ".fasta";
    public static final String CONTIG_MAPPINGS_FILE_ID = ".rna_contigs_mappings.tsv";

    // Junction anchors.
    // Minimum M flank (bp) for an N to count as a trusted splice junction. A flank below this carries too
    // little evidence to assert a junction (~1/4^8 chance a coincidental anchor matches), so an N with a
    // sub-threshold flank is treated as fabricated/untrustworthy: the discriminator won't swap the primary
    // onto it, and the rescue-fold / terminal-collapse passes re-evaluate the terminal anchor against the
    // genome rather than trusting bwa's placement. One threshold, used everywhere a junction anchor is judged.
    public static final int MIN_JUNCTION_ANCHOR = 8;
    public static final int ANNOTATED_JUNCTION_MIN_ANCHOR_BP = 1;
    public static final int ANNOTATED_JUNCTION_MIN_SOFTCLIP_ANCHOR_BP = 3;
    public static final int SPLICE_FLANKING_DELETION_MAX_BP = 5;

    // bwa-mem affine alignment scoring. The single source for every liftback pass that re-scores a boundary:
    // the BoundaryReclaim base walk (match/mismatch) and the LiftBackResolver CIGAR reconstruction (all four).
    // Soft-clips are not penalised (clipped bases score 0).
    public static final int MATCH = 1;
    public static final int MISMATCH = -4;
    public static final int GAP_OPEN = -6;
    public static final int GAP_EXTEND = -1;

    // MAPQ / AS thresholds.
    public static final int RESCUE_MAPQ = 60;
    public static final int SUPP_RESCUE_MAPQ_CAP = 55;
    public static final int PRIMARY_AS_UNMAP_THRESHOLD = 30;
    public static final int SUPP_AS_DROP_THRESHOLD = 30;

    // Junction canonicalize: max bp an intron is slid to reach a canonical motif. Bounds the search to the
    // size of a tx-contig boundary deletion (SPLICE_FLANKING_DELETION_MAX_BP), the only place it drifts.
    public static final int DEFAULT_MAX_SHIFT = 5;

    // Junction rescue defaults. Annotated-junction policy: min anchor overhang 3, min intron 21, max intron 1_000_000.
    public static final int DEFAULT_MIN_INTRON_LENGTH = 21;
    public static final int DEFAULT_MAX_INTRON_LENGTH = 1_000_000;
    public static final int DEFAULT_MIN_ANCHOR_OVERHANG = 3;
    public static final int DEFAULT_MAX_CHAIN_DEPTH = 4;
    // Overlap tolerance between primary and supp matched regions. BWA primary/supp often disagree
    // by 1-5 bases on junction placement; snapping to an annotated junction within this window recovers
    // the real splice. Snapping to the annotated junction is preferred; ties resolve via trust-primary.
    public static final int DEFAULT_SOFTCLIP_TOLERANCE = 5;
    // How many over-extended matched bases to trim when probing for an annotated donor/acceptor
    // (no-supp ref-verify path). BWA frequently extends past the true exon boundary into the intron;
    // measured over-extension caps at ~8 bases.
    public static final int DEFAULT_MAX_BOUNDARY_SHIFT = 8;
    // Minimum exon-proximal matched run for partial ref-verify rescue (softclip with an outer
    // unmatched tail). Partial matches need a longer anchor than full-match clips because the
    // unexplained residual raises false-positive risk; 11 makes coincidental matches implausible.
    public static final int DEFAULT_MIN_PARTIAL_MATCH_RUN = 11;

    // Tail-extension defaults. MaxExtension caps the walk to avoid consuming real unannotated junctions.
    public static final int DEFAULT_MIN_SOFTCLIP_LENGTH = 3;
    public static final int DEFAULT_MIN_EXTENSION = 3;
    public static final int DEFAULT_MAX_EXTENSION = 30;

    private TarsConstants() { }
}
