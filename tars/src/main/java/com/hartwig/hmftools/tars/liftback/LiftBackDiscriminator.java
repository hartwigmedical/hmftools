package com.hartwig.hmftools.tars.liftback;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

// Pure decision step on a record's full lifted alignment set (self + lifted XA alts).
// categorize() resolves the set to an Outcome (which source wins) plus a DecidingFeature (why); apply()
// mutates Dropped/IsPrimaryChoice for the features that need a swap (e.g. intronless ref alt vs. spliced
// parent gene). No I/O - separate from LiftBackResolver so the decision tree reads as one piece.
public final class LiftBackDiscriminator
{
    private LiftBackDiscriminator() { }

    // Per-alignment evidence emitted into TSV-A for post-hoc "why this outcome" analysis.
    public static class Features
    {
        public DecidingFeature Feature;
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
        Features features = new Features();

        if(alignments.isEmpty())
        {
            // Defensive: the primary path always supplies self, so this is unreachable in practice.
            features.Feature = null;
            return features;
        }

        Set<String> loci = new HashSet<>();
        Set<String> distinctCigars = new HashSet<>();
        boolean anyHasN = false;

        for(final LiftedAlignment alignment : alignments)
        {
            loci.add(locusKey(alignment));
            distinctCigars.add(alignment.LiftedCigar);
            if(alignment.cigarHasN())
            {
                anyHasN = true;
            }

            if(alignment.fromTxContig())
            {
                ++features.NumTxAlts;
                if(alignment.cigarHasRealNJunction())
                {
                    features.TxHasNCigar = true;
                }
                if(alignment.SoftClipAtBoundary)
                {
                    features.TxSoftClipAtBoundary = true;
                }
            }
            else
            {
                ++features.NumRefAlts;
                if(alignment.cigarHasSoftClip())
                {
                    features.RefSoftClipped = true;
                }
                else
                {
                    features.RefFullMatch = true;
                }
                if(alignment.cigarHasRealNJunction())
                {
                    features.RefHasNCigar = true;
                }
            }
        }

        boolean hasRef = features.NumRefAlts > 0;
        boolean hasTx = features.NumTxAlts > 0;

        if(loci.size() >= 2)
        {
            if(hasRef && hasTx)
            {
                // tx has an annotated junction; ref alts map intronlessly at other loci - usually an intronless
                // gene copy, though the rule does not verify that.
                if(features.TxHasNCigar && !features.RefHasNCigar)
                {
                    assign(features, DecidingFeature.JUNCTION_OVER_CONTIGUOUS);
                }
                else
                {
                    assign(features, DecidingFeature.MULTIMAPPER);
                }
            }
            else if(hasTx)
            {
                assign(features, DecidingFeature.SOLE_TX);
            }
            else
            {
                assign(features, DecidingFeature.SOLE_REF);
            }
            return features;
        }

        if(hasRef && !hasTx)
        {
            assign(features, DecidingFeature.SOLE_REF);
            return features;
        }
        if(hasTx && !hasRef)
        {
            assign(features, DecidingFeature.SOLE_TX);
            return features;
        }

        // single locus, both ref and tx - run discriminator
        if(distinctCigars.size() == 1 && !anyHasN)
        {
            assign(features, DecidingFeature.CONCORDANT);
        }
        else if(features.TxHasNCigar && features.RefSoftClipped)
        {
            assign(features, DecidingFeature.JUNCTION);
        }
        else if(features.TxHasNCigar && features.RefFullMatch && !features.RefSoftClipped)
        {
            // ref full-match across the supposed intron is strong evidence the read is unspliced
            // (pre-mRNA / retained intron / DNA contamination) - the tx N-CIGAR is the artifact.
            assign(features, DecidingFeature.REF_READS_THROUGH);
        }
        else if(features.TxSoftClipAtBoundary && features.RefFullMatch)
        {
            assign(features, DecidingFeature.REF_READS_THROUGH);
        }
        else
        {
            assign(features, DecidingFeature.AMBIGUOUS);
        }

        return features;
    }

    private static void assign(final Features features, final DecidingFeature feature)
    {
        features.Feature = feature;
    }

    public record ApplyResult(LiftedAlignment effectivePrimary, String note)
    {
    }

    // Mutates Dropped/IsPrimaryChoice and returns the effective primary and a short diagnostic note.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final DecidingFeature feature, final LiftedAlignment self)
    {
        if(feature == null)
        {
            return new ApplyResult(self, "");
        }
        return switch(feature)
        {
            case JUNCTION, JUNCTION_OVER_CONTIGUOUS -> promoteTxOverRef(alignments, self);
            case REF_READS_THROUGH -> promoteRefOverTx(alignments, self);
            default -> new ApplyResult(self, "");
        };
    }

    // Ref-favouring: if self is ref, drop tx alts. If self is tx, swap to the best ref alt.
    private static ApplyResult promoteRefOverTx(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(!self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(la.fromTxContig())
                {
                    la.Dropped = true;
                }
            return new ApplyResult(self, "");
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
        {
            return new ApplyResult(self, "");
        }

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la != self)
            {
                la.Dropped = true;
            }
        }
        return new ApplyResult(winner, "swapped_tx_to_ref");
    }

    // Tx-favouring: if self is tx, drop ref alts. If self is ref, swap to the best tx alt (self demoted to XA
    // to preserve the contiguous locus).
    private static ApplyResult promoteTxOverRef(final List<LiftedAlignment> alignments, final LiftedAlignment self)
    {
        if(self.fromTxContig())
        {
            for(final LiftedAlignment la : alignments)
                if(!la.fromTxContig())
                {
                    la.Dropped = true;
                }
            return new ApplyResult(self, "");
        }

        // Prefer the N-junction tx alt with fewest mismatches; fall back to any tx alt.
        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(la.fromTxContig() && la.cigarHasN())
            {
                if(winner == null || la.NumMismatches < winner.NumMismatches)
                {
                    winner = la;
                }
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
        {
            return new ApplyResult(self, "");
        }

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        for(final LiftedAlignment la : alignments)
        {
            if(!la.fromTxContig() && la != self)
            {
                la.Dropped = true;
            }
        }
        return new ApplyResult(winner, "swapped_ref_to_tx");
    }

    private static String locusKey(final LiftedAlignment la)
    {
        return la.LiftedChrom + ":" + la.LiftedPos;
    }
}
