package com.hartwig.hmftools.redux.splice;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

// pure decision step on a record's full lifted alignment set (self + lifted XA alts):
//   1. categorize() walks the alignments and returns the LiftBackCategory along with the per-alignment
//      evidence the resolver feeds into LiftBackResult / TSV-A.
//   2. apply() resolves the discriminator outcome for the categories that need one: marks losing-side
//      alts Dropped, and may swap the BAM primary onto a winning XA alt (e.g. when bwa picked an
//      intronless paralog over the spliced parent gene).
//
// No SAMRecord, no I/O — kept separate from LiftBackResolver so the decision tree reads as one piece.
public final class LiftBackDiscriminator
{
    // minimum flanking M length for an N operator to be treated as a real splice junction in
    // discriminator decisions. Lift-back of a tx alignment whose last bases bleed into the next exon
    // can produce N with anchors of 1-3 bp — those are not evidence of splicing and must not swap the
    // primary off a clean ref full-match. 8 bp gives ~1-in-65k random-match probability.
    static final int MIN_REAL_N_ANCHOR = 8;

    private LiftBackDiscriminator() {}

    // per-alignment evidence the discriminator collected, emitted into TSV-A so post-hoc analysis can
    // answer "why did this read land in category X?".
    public static class Features
    {
        public LiftBackCategory Category;
        public int NumRefAlts;
        public int NumTxAlts;
        public boolean TxHasNCigar;
        public boolean TxSoftClipAtBoundary;
        public boolean RefSoftClipped;
        public boolean RefFullMatch;
        public boolean RefHasNCigar;
    }

    public static Features categorize(final List<LiftedAlignment> alignments)
    {
        final Features features = new Features();

        if(alignments.isEmpty())
        {
            features.Category = LiftBackCategory.LIFT_FAILED;
            return features;
        }

        final Set<String> loci = new HashSet<>();
        final Set<String> distinctCigars = new HashSet<>();
        boolean anyHasN = false;

        for(final LiftedAlignment alignment : alignments)
        {
            loci.add(locusKey(alignment));
            distinctCigars.add(alignment.LiftedCigar);
            if(alignment.cigarHasN())
                anyHasN = true;

            if(alignment.fromTxContig())
            {
                ++features.NumTxAlts;
                if(alignment.cigarHasRealNJunction(MIN_REAL_N_ANCHOR))
                    features.TxHasNCigar = true;
                if(alignment.SoftClipAtBoundary)
                    features.TxSoftClipAtBoundary = true;
            }
            else
            {
                ++features.NumRefAlts;
                if(alignment.cigarHasSoftClip())
                    features.RefSoftClipped = true;
                else
                    features.RefFullMatch = true;
                if(alignment.cigarHasRealNJunction(MIN_REAL_N_ANCHOR))
                    features.RefHasNCigar = true;
            }
        }

        final boolean hasRef = features.NumRefAlts > 0;
        final boolean hasTx = features.NumTxAlts > 0;

        if(loci.size() >= 2)
        {
            if(hasRef && hasTx)
            {
                // tx found an annotated junction at one locus, ref alts hit other loci intronlessly —
                // almost always processed pseudogenes / paralogs. Favour tx and swap the BAM primary.
                if(features.TxHasNCigar && !features.RefHasNCigar)
                    features.Category = LiftBackCategory.BOTH_MULTI_TX_JUNCTION;
                else
                    features.Category = LiftBackCategory.BOTH_MULTI;
            }
            else if(hasTx)
                features.Category = LiftBackCategory.TX_MULTI;
            else
                features.Category = LiftBackCategory.REF_MULTI;
            return features;
        }

        // single locus
        if(hasRef && !hasTx)
        {
            features.Category = LiftBackCategory.REF_SINGLE;
            return features;
        }
        if(hasTx && !hasRef)
        {
            features.Category = LiftBackCategory.TX_SINGLE;
            return features;
        }

        // single locus, both ref and tx contributed — run discriminator
        if(distinctCigars.size() == 1 && !anyHasN)
            features.Category = LiftBackCategory.BOTH_AGREE;
        else if(features.TxHasNCigar && features.RefSoftClipped)
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP;
        else if(features.TxHasNCigar && features.RefFullMatch && !features.RefSoftClipped)
            // ref full-match across the supposed intron with NM=0 is overwhelming evidence the read is
            // genuinely unspliced (pre-mRNA / retained-intron / DNA contamination). The tx N-CIGAR is
            // the artifact — favour ref.
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH;
        else if(features.TxSoftClipAtBoundary && features.RefFullMatch)
            features.Category = LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH;
        else
            features.Category = LiftBackCategory.BOTH_AMBIGUOUS;

        return features;
    }

    public record Outcome(LiftedAlignment effectivePrimary, String note) {}

    // applies confident discriminator outcomes by mutating LiftedAlignment.Dropped / IsPrimaryChoice.
    // Returns the effective primary (the alignment whose lifted coords become the BAM record's coords)
    // and a short note for TSV diagnostics.
    public static Outcome apply(
            final List<LiftedAlignment> alignments, final LiftBackCategory category, final LiftedAlignment self)
    {
        switch(category)
        {
            case BOTH_TX_JUNCTION_REF_SOFTCLIP:
            case BOTH_MULTI_TX_JUNCTION:
                return promoteTxOverRef(alignments, self);

            case BOTH_TX_JUNCTION_REF_MATCH:
                return promoteRefOverTx(alignments, self);

            case BOTH_TX_SOFTCLIP_REF_MATCH:
                if(!self.fromTxContig())
                {
                    for(final LiftedAlignment la : alignments)
                        if(la.fromTxContig())
                            la.Dropped = true;
                    return new Outcome(self, "");
                }
                return new Outcome(self, "self_was_tx_no_swap");

            default:
                return new Outcome(self, "");
        }
    }

    // ref-favouring outcome. If self is ref, drop all tx alts. If self is tx, swap: promote the winning
    // ref alt to primary, demote self (tx) to XA, drop other tx alts.
    private static Outcome promoteRefOverTx(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(!self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(la.fromTxContig())
                    la.Dropped = true;
            return new Outcome(self, "");
        }

        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(!la.fromTxContig())
            {
                winner = la;
                break;
            }
        }
        if(winner == null)
            return new Outcome(self, "");

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la != self)
                la.Dropped = true;
        }
        return new Outcome(winner, "swapped_tx_to_ref");
    }

    // tx-favouring outcome. If self is tx, drop all ref alts. If self is ref, swap: promote the winning
    // tx alt to primary, demote self to XA (preserves the original locus as an informative paralog hit),
    // and drop other ref alts.
    private static Outcome promoteTxOverRef(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(!la.fromTxContig())
                    la.Dropped = true;
            return new Outcome(self, "");
        }

        // prefer the tx alt that actually carries the N junction; fall back to any tx alt
        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la.cigarHasN())
            {
                winner = la;
                break;
            }
        }
        if(winner == null)
        {
            for(final LiftedAlignment la : alignments)
            {
                if(la.fromTxContig())
                {
                    winner = la;
                    break;
                }
            }
        }
        if(winner == null)
            return new Outcome(self, "");

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        for(final LiftedAlignment la : alignments)
        {
            if(!la.fromTxContig() && la != self)
                la.Dropped = true;
        }
        return new Outcome(winner, "swapped_ref_to_tx");
    }

    private static String locusKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos;
    }
}
