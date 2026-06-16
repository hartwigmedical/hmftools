package com.hartwig.hmftools.tars.liftback.rescue;

// Drives per-decision counters in RescueStatistics.
public enum RescueRejectReason
{
    NO_TERMINAL_SOFTCLIP,         // primary cigar has no leading or trailing S to extend across
    NO_MATCHING_SUPP,             // no supplementary had a complementary cigar shape
    DIFFERENT_CHROMOSOME,
    OPPOSITE_STRAND,
    READ_COVERAGE_OVERLAP,        // primary's M + supp's M overlap on the read
    READ_COVERAGE_GAP,            // primary's M + supp's M leave a gap on the read
    INTRON_TOO_SHORT,
    INTRON_TOO_LONG,
    SHORT_ANCHOR,                 // primary or supp matched portion < MinAnchorOverhang
    NOVEL_JUNCTION,               // candidate intron not in annotated set (and AnnotatedOnly is true)
    COMPLEX_CIGAR_SHAPE,          // hard clip, indel adjacent to softclip boundary, etc.
    AMBIGUOUS_SUPP_CHOICE,        // multiple equally-good candidate supps; refuse to guess
    LOW_PRIMARY_MAPQ,             // primary MAPQ below the merge floor - placement not trusted enough
    MULTIPLE_SUPPS_IN_REACH,      // more than one supp within merge reach; refuse to guess which splice
    // ref-verify path: primary has terminal softclip but no matching supp
    REF_VERIFY_NO_CANDIDATE_EXON, // no annotated junction's adjacent exon lines up with the softclip
    REF_VERIFY_MISMATCH_TOO_HIGH, // softclipped read bases don't match the candidate exon's ref bases
    REF_VERIFY_NO_REF_SOURCE,     // ref-verify requested but no RefSequenceSource configured
    REF_VERIFY_AMBIGUOUS,         // multiple annotated downstream/upstream exons match; refuse to guess
    REF_VERIFY_SHORT_PARTIAL_RUN  // partial match (outer residual left clipped) shorter than MinPartialMatchRun
}
