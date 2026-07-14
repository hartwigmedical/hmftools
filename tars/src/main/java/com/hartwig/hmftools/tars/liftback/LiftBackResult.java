package com.hartwig.hmftools.tars.liftback;

import java.util.List;

// Per-record decision result. Drives the BAM record rewrite and one TSV-A row.
// recordState says whether the record reached the discriminator; for RESOLVED records outcome/decidingFeature
// carry the ref-vs-tx decision and its reason. liftedAlignments holds the full alignment set (self + lifted
// XA alts) for TSV-B.
public record LiftBackResult(
        RecordState recordState,
        DecidingFeature decidingFeature,
        boolean swapped,
        Composition comp,
        RecordRole role,
        String finalChrom,
        int finalPos,
        String finalCigar,
        boolean negativeStrand,
        boolean hasNCigar,
        int inputMapq,
        int updatedMapq,
        int numXaAlts,
        int numRefAlts,
        int numTxAlts,
        int numLoci,
        int numDistinctCigarsAtPrimaryLocus,
        boolean txHasNCigar,
        boolean txSoftClipAtBoundary,
        boolean refSoftClipped,
        boolean refFullMatch,
        String geneIds,
        String notes,
        // +1/-1 for tx-contig-derived primaries; 0 for ref-only or supplementary-resolve/tail-extend origins.
        // Used by the writer to set XS:A:+/- on spliced records.
        int transcriptStrand,
        List<LiftedAlignment> liftedAlignments)
{
    // Outcome is derived from the deciding feature rather than stored alongside it; non-decision records
    // (null feature: unmapped / lift-failed / supplementary) are UNRESOLVED.
    public Outcome outcome()
    {
        return decidingFeature != null ? decidingFeature.outcome() : Outcome.UNRESOLVED;
    }

    // Copy with a post-lift cigar pass applied: position, cigar, N-flag, MAPQ and an appended note change;
    // every other field is preserved. Shared by the supplementary-resolve / collapse / tail-extend / canonicalize passes.
    public LiftBackResult withRevisedCigar(
            final int newPos, final String newCigar, final boolean newHasNCigar, final int newUpdatedMapq, final String note)
    {
        return new LiftBackResult(
                recordState, decidingFeature, swapped, comp, role,
                finalChrom, newPos, newCigar,
                negativeStrand, newHasNCigar,
                inputMapq, newUpdatedMapq,
                numXaAlts, numRefAlts, numTxAlts,
                numLoci, numDistinctCigarsAtPrimaryLocus,
                txHasNCigar, txSoftClipAtBoundary,
                refSoftClipped, refFullMatch,
                geneIds,
                appendNote(notes, note),
                transcriptStrand,
                liftedAlignments);
    }

    // UNMAPPED result with the given role and a diagnostic note. Used when a lifted placement is rejected
    // post-lift (e.g. it falls in an excluded region) and the record must be flipped to unmapped.
    public static LiftBackResult unmapped(final RecordRole role, final String note)
    {
        return new LiftBackResult(
                RecordState.UNMAPPED, null, false, Composition.NONE, role,
                "*", 0, "*",
                false, false, 0, 0,
                0, 0, 0,
                0, 0,
                false, false, false, false,
                "", note,
                0,
                List.of());
    }

    private static String appendNote(final String existing, final String note)
    {
        if(existing == null || existing.isEmpty())
        {
            return note;
        }
        return existing + ";" + note;
    }

    public enum RecordRole
    {
        PRIMARY,
        SUPPLEMENTARY
    }

    // post-drop view drives TSV-A + XA rebuild; pre-drop view drives LiftBackStats summary
    public enum Composition
    {
        REF_ONLY,
        TX_ONLY,
        REF_AND_TX,
        NONE;

        public static Composition fromAlignments(final List<LiftedAlignment> alignments)
        {
            if(alignments.isEmpty())
            {
                return NONE;
            }

            boolean hasRef = false;
            boolean hasTx = false;
            for(LiftedAlignment alignment : alignments)
            {
                if(alignment.fromTxContig())
                {
                    hasTx = true;
                }
                else
                {
                    hasRef = true;
                }
            }

            if(hasRef && hasTx)
            {
                return REF_AND_TX;
            }
            if(hasTx)
            {
                return TX_ONLY;
            }
            return REF_ONLY;
        }
    }
}
