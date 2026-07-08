package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
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
            // an alt the overhang gate collapsed to a contiguous alignment is marked Dropped before the
            // discriminator runs; it is a fabricated placement, so it must not count toward the features.
            if(alignment.Dropped)
            {
                continue;
            }

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

    // Deterministic overload: respects bwa's order (no randomization). Used where the caller has no read seed.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final DecidingFeature feature, final LiftedAlignment self)
    {
        return apply(alignments, feature, self, 0, true);
    }

    // Mutates Dropped/IsPrimaryChoice and returns the effective primary and a short diagnostic note.
    // When bwa expressed no priority (bwaHasPriority false, ie MAPQ 0) the not-confident outcomes first rank the
    // candidates by their recomputed genomic score (set pre-discriminator) and take the clear winner: the
    // highest-scoring locus among multi-mappers, or the higher-scoring side of an ambiguous tx-vs-ref call. Only
    // when the top candidates tie on score - OR the candidates were never scored (e.g. a split read left for
    // supplementary-resolve) - does it fall back to bwa's seed-based pseudo-random pick. bwaHasPriority true leaves
    // bwa's order untouched.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final DecidingFeature feature, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority)
    {
        if(feature == null)
        {
            return new ApplyResult(self, "");
        }
        return switch(feature)
        {
            case JUNCTION, JUNCTION_OVER_CONTIGUOUS -> promoteTxOverRef(alignments, self);
            case REF_READS_THROUGH -> promoteRefOverTx(alignments, self);
            case AMBIGUOUS -> bwaHasPriority ? new ApplyResult(self, "") : randomTxVsRef(alignments, self, seed);
            case MULTIMAPPER, SOLE_TX, SOLE_REF -> bwaHasPriority ? new ApplyResult(self, "") : randomLocus(alignments, self, seed);
            default -> new ApplyResult(self, "");
        };
    }

    // Pick the primary locus among multi-mappers: the highest-scoring locus, with the seeded coin only among loci
    // tied at the top. When no candidate was scored (a split read whose scoring is skipped so supplementary-resolve
    // can reconstruct it), fall back to bwa's exact behaviour - a seed pick over all loci, first placement there -
    // so the merge path is untouched. Every non-dropped placement stays (rides in XA); single-locus sets keep self.
    private static ApplyResult randomLocus(final List<LiftedAlignment> alignments, final LiftedAlignment self, final int seed)
    {
        List<String> loci = new ArrayList<>();
        for(final LiftedAlignment la : alignments)
        {
            if(la.Dropped)
            {
                continue;
            }
            String key = locusKey(la);
            if(!loci.contains(key))
            {
                loci.add(key);
            }
        }
        if(loci.size() < 2)
        {
            return new ApplyResult(self, "");
        }

        int topScore = Integer.MIN_VALUE;
        for(final String locus : loci)
        {
            topScore = Math.max(topScore, bestScoreAtLocus(alignments, locus));
        }
        boolean scored = topScore != Integer.MIN_VALUE;

        List<String> pickable = new ArrayList<>();
        for(final String locus : loci)
        {
            if(!scored || bestScoreAtLocus(alignments, locus) == topScore)
            {
                pickable.add(locus);
            }
        }

        boolean tie = pickable.size() > 1;
        String pickedLocus = pickable.get(Math.floorMod(seed, pickable.size()));
        if(pickedLocus.equals(locusKey(self)))
        {
            return new ApplyResult(self, "");
        }

        // unscored: first placement at the locus (bwa behaviour); scored: the best-scoring placement there.
        LiftedAlignment winner = null;
        for(final LiftedAlignment la : alignments)
        {
            if(la.Dropped || !locusKey(la).equals(pickedLocus))
            {
                continue;
            }
            if(winner == null || (scored && outranks(la, winner)))
            {
                winner = la;
            }
        }
        if(winner == null)
        {
            return new ApplyResult(self, "");
        }

        self.IsPrimaryChoice = false;
        winner.IsPrimaryChoice = true;
        return new ApplyResult(winner, (scored && !tie) ? "score_locus" : "random_locus");
    }

    private static int bestScoreAtLocus(final List<LiftedAlignment> alignments, final String locus)
    {
        int best = Integer.MIN_VALUE;
        for(final LiftedAlignment la : alignments)
        {
            if(!la.Dropped && locusKey(la).equals(locus))
            {
                best = Math.max(best, la.GenomicScore);
            }
        }
        return best;
    }

    // Single locus with both a tx and a ref placement and no shape rule to separate them: keep the higher genomic
    // score, or a seeded coin when they tie (or were not scored). The loser is a different cigar at the same locus,
    // kept in XA (not dropped).
    private static ApplyResult randomTxVsRef(final List<LiftedAlignment> alignments, final LiftedAlignment self, final int seed)
    {
        LiftedAlignment tx = bestOfSource(alignments, true);
        LiftedAlignment ref = bestOfSource(alignments, false);
        if(tx == null || ref == null)
        {
            return new ApplyResult(self, "");
        }

        LiftedAlignment winner;
        String note;
        if(bothScored(tx, ref) && tx.GenomicScore != ref.GenomicScore)
        {
            // a real genomic score separates them, so it is not a true tie: keep the higher-scoring side.
            winner = tx.GenomicScore > ref.GenomicScore ? tx : ref;
            note = winner == tx ? "score_tx" : "score_ref";
        }
        else
        {
            // genuinely tied (or scores unavailable): fall back to the read-seeded coin so ties stay reproducible.
            boolean pickTx = (seed & 1) == 0;
            winner = pickTx ? tx : ref;
            note = pickTx ? "random_tx" : "random_ref";
        }
        for(final LiftedAlignment la : alignments)
        {
            la.IsPrimaryChoice = la == winner;
        }
        return new ApplyResult(winner, note);
    }

    // Best non-dropped alignment from the tx side (fromTx true) or ref side (fromTx false).
    private static LiftedAlignment bestOfSource(final List<LiftedAlignment> alignments, final boolean fromTx)
    {
        LiftedAlignment best = null;
        for(final LiftedAlignment la : alignments)
        {
            if(la.Dropped || la.fromTxContig() != fromTx)
            {
                continue;
            }
            if(best == null || outranks(la, best))
            {
                best = la;
            }
        }
        return best;
    }

    // Higher recomputed genomic score wins; an equal score (or scores not computed) falls back to fewer mismatches.
    private static boolean outranks(final LiftedAlignment a, final LiftedAlignment b)
    {
        if(a.GenomicScore != b.GenomicScore)
        {
            return a.GenomicScore > b.GenomicScore;
        }
        return a.NumMismatches < b.NumMismatches;
    }

    private static boolean bothScored(final LiftedAlignment a, final LiftedAlignment b)
    {
        return a.GenomicScore != Integer.MIN_VALUE && b.GenomicScore != Integer.MIN_VALUE;
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
