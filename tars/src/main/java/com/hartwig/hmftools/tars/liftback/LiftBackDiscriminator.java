package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.MATE_PROXIMITY_MAX_DISTANCE;

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

        for(LiftedAlignment alignment : alignments)
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

    // Mate-agnostic overload: single-end reads and callers with no lifted mate.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final DecidingFeature feature, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority)
    {
        return apply(alignments, feature, self, seed, bwaHasPriority, null);
    }

    // Mutates Dropped/IsPrimaryChoice and returns the effective primary and a short diagnostic note.
    // When bwa expressed no priority (bwaHasPriority false, ie MAPQ 0) the not-confident outcomes first rank the
    // candidates by their recomputed genomic score (set pre-discriminator) and take the clear winner: the
    // highest-scoring locus among multi-mappers, or the higher-scoring side of an ambiguous tx-vs-ref call. Only
    // when the top candidates tie on score - OR the candidates were never scored (e.g. a split read left for
    // supplementary-resolve) - does it fall back to the mate-proximity / junction / seed tie-breaks. bwaHasPriority
    // true leaves bwa's order untouched. mateInfo is the partner mate's lifted primary placement, or null.
    public static ApplyResult apply(
            final List<LiftedAlignment> alignments, final DecidingFeature feature, final LiftedAlignment self,
            final int seed, final boolean bwaHasPriority, final LiftedMateInfo mateInfo)
    {
        if(feature == null || feature == DecidingFeature.CONCORDANT || bwaHasPriority)
        {
            return new ApplyResult(self, "");
        }
        return pickByScore(alignments, self, seed, mateInfo);
    }

    // Highest recomputed genome score wins ("score"). Top-score ties are settled in order by: mate proximity ("mate"),
    // then a spliced placement over a same-locus soft-clip ("junction"), then a read-name seed ("random"). Nothing is
    // dropped; losers ride in XA. Unscored candidates (a split read left for Step 3) keep bwa's primary.
    private static ApplyResult pickByScore(
            final List<LiftedAlignment> alignments, final LiftedAlignment self, final int seed, final LiftedMateInfo mateInfo)
    {
        List<LiftedAlignment> candidates = new ArrayList<>();
        for(LiftedAlignment alignment : alignments)
        {
            if(!alignment.Dropped)
            {
                candidates.add(alignment);
            }
        }
        if(candidates.size() < 2)
        {
            return new ApplyResult(self, "");
        }

        int topScore = Integer.MIN_VALUE;
        for(LiftedAlignment alignment : candidates)
        {
            topScore = Math.max(topScore, alignment.GenomicScore);
        }
        if(topScore == Integer.MIN_VALUE)
        {
            return new ApplyResult(self, ""); // unscored (split read left for Step 3): keep bwa's placement
        }

        // Collapse identical placements (same locus + CIGAR from different sources, e.g. a ref self and a tx alt that
        // lift to the same contiguous alignment) so the tie is over distinct placements, not weighted by source count.
        List<LiftedAlignment> top = new ArrayList<>();
        Set<String> topKeys = new HashSet<>();
        for(LiftedAlignment alignment : candidates)
        {
            if(alignment.GenomicScore == topScore && topKeys.add(placementKey(alignment)))
            {
                top.add(alignment);
            }
        }

        boolean tie = top.size() > 1;
        LiftedAlignment winner;
        String note;
        if(!tie)
        {
            winner = top.get(0);
            note = "score";
        }
        else
        {
            List<LiftedAlignment> contenders = mateProximalSubset(top, mateInfo);
            if(contenders.size() == 1)
            {
                winner = contenders.get(0);
                note = "mate";
            }
            else
            {
                LiftedAlignment junction = preferJunctionOverSoftClip(contenders);
                if(junction != null)
                {
                    winner = junction;
                    note = "junction";
                }
                else
                {
                    winner = contenders.get(Math.floorMod(seed, contenders.size()));
                    note = "random";
                }
            }
        }
        for(LiftedAlignment alignment : alignments)
        {
            alignment.IsPrimaryChoice = alignment == winner;
        }
        return new ApplyResult(winner, note);
    }

    // Tied candidates proximal to the mate; the full set unchanged when the mate is absent or does not discriminate.
    private static List<LiftedAlignment> mateProximalSubset(final List<LiftedAlignment> top, final LiftedMateInfo mateInfo)
    {
        if(mateInfo == null || mateInfo.unmapped() || mateInfo.chromosome() == null)
        {
            return top;
        }
        List<LiftedAlignment> near = new ArrayList<>();
        for(LiftedAlignment alignment : top)
        {
            if(isMateProximal(alignment, mateInfo))
            {
                near.add(alignment);
            }
        }
        return (near.isEmpty() || near.size() == top.size()) ? top : near;
    }

    private static boolean isMateProximal(final LiftedAlignment alignment, final LiftedMateInfo mateInfo)
    {
        if(!mateInfo.chromosome().equals(alignment.LiftedChrom))
        {
            return false;
        }
        int gap;
        if(alignment.LiftedPos > mateInfo.alignmentEnd())
        {
            gap = alignment.LiftedPos - mateInfo.alignmentEnd();
        }
        else if(alignment.LiftedPos < mateInfo.alignmentStart())
        {
            gap = mateInfo.alignmentStart() - alignment.LiftedPos;
        }
        else
        {
            gap = 0;
        }
        return gap <= MATE_PROXIMITY_MAX_DISTANCE;
    }

    // Tie-break within an equal-top-score set: a spliced placement (a real N junction) beats a clipped placement
    // (soft-clip, no N) at the same lifted locus. Both describe the same read at the same start; the junction is the
    // correct RNA interpretation - bwa soft-clipped rather than cross the intron - so it is not left to the coin.
    private static LiftedAlignment preferJunctionOverSoftClip(final List<LiftedAlignment> top)
    {
        for(LiftedAlignment junction : top)
        {
            if(!junction.cigarHasRealNJunction())
            {
                continue;
            }
            for(LiftedAlignment clipped : top)
            {
                if(clipped != junction
                        && locusKey(clipped).equals(locusKey(junction))
                        && clipped.cigarHasSoftClip()
                        && !clipped.cigarHasRealNJunction())
                {
                    return junction;
                }
            }
        }
        return null;
    }

    private static String locusKey(final LiftedAlignment alignment)
    {
        return alignment.LiftedChrom + ":" + alignment.LiftedPos;
    }

    private static String placementKey(final LiftedAlignment alignment)
    {
        return locusKey(alignment) + ":" + alignment.LiftedCigar;
    }
}
