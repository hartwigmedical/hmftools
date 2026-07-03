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
    // onto it (LiftedAlignment.cigarHasRealNJunction). The overhang gate scores terminal anchors against the
    // genome instead of applying this length threshold.
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

    // Overhang gate: keep a terminal splice-junction anchor only if its own affine alignment score (matches +1,
    // mismatches -4) at the far exon exceeds this. Equals the bwa soft-clip penalty, so a junction is kept only
    // when aligning the anchor scores better than clipping it (a clean 6bp anchor survives, a 5bp anchor does not).
    public static final int MIN_OVERHANG_SCORE = 5;

    // Overhang gate length trust: an overhang longer than this is trusted outright and never scored/collapsed - a long
    // anchor is a real junction even with a couple of mismatches (SNVs). Only short overhangs (<= this) run the score
    // check above. Applies to the terminal-anchor collapse in both the single-junction and multi-junction cases.
    public static final int MIN_OVERHANG_LENGTH = 12;

    // MAPQ / AS thresholds.
    // Confident placement MAPQ: bumped onto a single-locus MAPQ-0 primary (decidePrimaryMapq) and onto a
    // supplementary-resolve merge whose components were all MAPQ 0 (the merge is the placement evidence).
    public static final int CONFIDENT_MAPQ = 60;
    public static final int PRIMARY_AS_UNMAP_THRESHOLD = 30;
    public static final int SUPP_AS_DROP_THRESHOLD = 30;

    // Supplementary-resolve defaults. Annotated-junction policy: min intron 21, max intron 1_000_000.
    public static final int DEFAULT_MIN_INTRON_LENGTH = 21;
    public static final int DEFAULT_MAX_INTRON_LENGTH = 1_000_000;

    // max supplementaries merged into one read: 1 (typical) or 2 (a middle exon clipped both ends). Deeper chains
    // are vanishingly rare (depth 3-4 was 72 reads in a whole sample), so 2 covers essentially all real cases.
    public static final int DEFAULT_MAX_SUPP_MERGES = 2;

    // Overlap tolerance (bp) between primary and supp matched regions in a supplementary-resolve merge. BWA primary/supp often
    // disagree by 1-5 bases on junction placement; the snap resolves within this window. Compile-time only (not a flag).
    public static final int DEFAULT_SOFTCLIP_TOLERANCE = 5;

    // How many over-extended matched bases to trim when probing for an annotated donor/acceptor
    // (no-supp ref-verify path). BWA frequently extends past the true exon boundary into the intron;
    // measured over-extension caps at ~8 bases.
    public static final int DEFAULT_MAX_BOUNDARY_SHIFT = 8;

    private TarsConstants() { }
}
