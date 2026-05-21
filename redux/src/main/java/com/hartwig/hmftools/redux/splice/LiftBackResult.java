package com.hartwig.hmftools.redux.splice;

import java.util.List;

// per-record decision result. Drives the BAM record rewrite and one TSV-A row.
// liftedAlignments contains every component of the record's alignment set (self + lifted XA alts) for TSV-B.
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
        List<LiftedAlignment> liftedAlignments)
{
    public enum RecordRole
    {
        PRIMARY,
        SUPPLEMENTARY
    }

    // contig composition of an alignment set. Two views: (a) post-drop "kept" alignments drives TSV-A +
    // XA rebuild; (b) pre-drop full alignment set drives the LiftBackStats summary.
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
