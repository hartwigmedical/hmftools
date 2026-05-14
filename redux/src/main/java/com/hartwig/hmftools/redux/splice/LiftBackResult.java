package com.hartwig.hmftools.redux.splice;

import java.util.List;

// per-record decision result. Drives the BAM record rewrite and one TSV-A row.
// LiftedAlignments contains every component of the record's alignment set (self + lifted XA alts) for TSV-B.
public class LiftBackResult
{
    public enum RecordRole
    {
        PRIMARY,
        SUPPLEMENTARY
    }

    // describes the contig composition of an alignment set (the BAM record's self + lifted XA alts).
    // used in two views: (a) post-drop "kept" alignments (LiftBackResult.Comp) drives TSV-A + XA rebuild;
    // (b) pre-drop full alignment set, derived on the fly from LiftedAlignments, drives the LiftBackStats summary.
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
            for(LiftedAlignment alignment : alignments)
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

    public final LiftBackCategory Category;
    public final Composition Comp;
    public final RecordRole Role;
    public final String FinalChrom;
    public final int FinalPos;
    public final String FinalCigar;
    public final boolean NegativeStrand;
    public final boolean HasNCigar;
    public final int BwaMapq;
    public final int UpdatedMapq;
    public final int NumXaAlts;
    public final int NumRefAlts;
    public final int NumTxAlts;
    public final int NumLoci;
    public final int NumDistinctCigarsAtPrimaryLocus;
    public final boolean TxHasNCigar;
    public final boolean TxSoftClipAtBoundary;
    public final boolean RefSoftClipped;
    public final boolean RefFullMatch;
    public final String GeneIds;
    public final String Notes;
    public final List<LiftedAlignment> LiftedAlignments;

    public LiftBackResult(
            final LiftBackCategory category, final Composition composition, final RecordRole role,
            final String finalChrom, final int finalPos, final String finalCigar, final boolean negativeStrand,
            final boolean hasNCigar, final int bwaMapq, final int updatedMapq,
            final int numXaAlts, final int numRefAlts, final int numTxAlts,
            final int numLoci, final int numDistinctCigarsAtPrimaryLocus,
            final boolean txHasNCigar, final boolean txSoftClipAtBoundary,
            final boolean refSoftClipped, final boolean refFullMatch,
            final String geneIds, final String notes,
            final List<LiftedAlignment> liftedAlignments)
    {
        Category = category;
        Comp = composition;
        Role = role;
        FinalChrom = finalChrom;
        FinalPos = finalPos;
        FinalCigar = finalCigar;
        NegativeStrand = negativeStrand;
        HasNCigar = hasNCigar;
        BwaMapq = bwaMapq;
        UpdatedMapq = updatedMapq;
        NumXaAlts = numXaAlts;
        NumRefAlts = numRefAlts;
        NumTxAlts = numTxAlts;
        NumLoci = numLoci;
        NumDistinctCigarsAtPrimaryLocus = numDistinctCigarsAtPrimaryLocus;
        TxHasNCigar = txHasNCigar;
        TxSoftClipAtBoundary = txSoftClipAtBoundary;
        RefSoftClipped = refSoftClipped;
        RefFullMatch = refFullMatch;
        GeneIds = geneIds;
        Notes = notes;
        LiftedAlignments = liftedAlignments;
    }
}
