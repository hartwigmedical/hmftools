package com.hartwig.hmftools.redux.splice;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

// Pure decision step on a record's full lifted alignment set (self + lifted XA alts).
// categorize() classifies the set; apply() mutates Dropped/IsPrimaryChoice for categories that need
// a swap (e.g. intronless paralog vs. spliced parent gene). No I/O — separate from LiftBackResolver
// so the decision tree reads as one piece.
public final class LiftBackDiscriminator
{
    // Minimum flanking M length for an N to count as a real splice junction. Short anchors (1-3 bp)
    // from tx-boundary bleed are not splicing evidence; 8 bp gives ~1-in-65k random-match probability.
    static final int MIN_REAL_N_ANCHOR = 8;

    private LiftBackDiscriminator() {}

    // Per-alignment evidence emitted into TSV-A for post-hoc "why this category" analysis.
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
                // tx has annotated junction, ref alts map intronlessly at other loci — typically paralogs/pseudogenes.
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

        // single locus, both ref and tx — run discriminator
        if(distinctCigars.size() == 1 && !anyHasN)
            features.Category = LiftBackCategory.BOTH_AGREE;
        else if(features.TxHasNCigar && features.RefSoftClipped)
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_SOFTCLIP;
        else if(features.TxHasNCigar && features.RefFullMatch && !features.RefSoftClipped)
            // ref full-match across the supposed intron is strong evidence the read is unspliced
            // (pre-mRNA / retained intron / DNA contamination) — the tx N-CIGAR is the artifact.
            features.Category = LiftBackCategory.BOTH_TX_JUNCTION_REF_MATCH;
        else if(features.TxSoftClipAtBoundary && features.RefFullMatch)
            features.Category = LiftBackCategory.BOTH_TX_SOFTCLIP_REF_MATCH;
        else
            features.Category = LiftBackCategory.BOTH_AMBIGUOUS;

        return features;
    }

    public record Outcome(LiftedAlignment effectivePrimary, String note) {}

    // Mutates Dropped/IsPrimaryChoice and returns the effective primary and a short diagnostic note.
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

    // Ref-favouring: if self is ref, drop tx alts. If self is tx, swap to the best ref alt.
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

    // Tx-favouring: if self is tx, drop ref alts. If self is ref, swap to the best tx alt (self demoted to XA to preserve the paralog locus).
    private static Outcome promoteTxOverRef(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(!la.fromTxContig())
                    la.Dropped = true;
            return new Outcome(self, "");
        }

        // Prefer the N-junction tx alt with fewest mismatches; fall back to any tx alt.
        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la.cigarHasN())
            {
                if(winner == null || la.NumMismatches < winner.NumMismatches)
                    winner = la;
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
