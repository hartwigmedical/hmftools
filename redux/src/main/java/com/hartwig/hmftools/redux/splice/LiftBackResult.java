package com.hartwig.hmftools.redux.splice;

import java.util.List;

// Per-record decision result. Drives the BAM record rewrite and one TSV-A row.
// liftedAlignments holds the full alignment set (self + lifted XA alts) for TSV-B.
public record LiftBackResult(
        LiftBackCategory category,
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
        // +1/−1 for tx-contig-derived primaries; 0 for ref-only or rescue/tail-extend origins.
        // Used by the writer to set XS:A:+/- on spliced records.
        int transcriptStrand,
        List<LiftedAlignment> liftedAlignments)
{
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
                return NONE;

            boolean hasRef = false;
            boolean hasTx = false;
            for(final LiftedAlignment alignment : alignments)
            {
                if(alignment.fromTxContig())
                    hasTx = true;
                else
                    hasRef = true;
            }

            if(hasRef && hasTx)
                return REF_AND_TX;
            if(hasTx)
                return TX_ONLY;
            return REF_ONLY;
        }
    }
}
