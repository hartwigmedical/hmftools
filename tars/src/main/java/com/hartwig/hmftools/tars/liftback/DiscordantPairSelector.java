package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.tars.liftback.features.PlacementPairPriority;

// Selects the two emitted placements as one fragment when either mate has MAPQ 0.
final class DiscordantPairSelector
{
    private DiscordantPairSelector() { }

    static Optional<Choice> select(final MatePlacements first, final MatePlacements second, final int seed)
    {
        List<PairCandidate> bestPairs = new ArrayList<>();
        PlacementPairPriority bestPriority = null;

        for(LiftedAlignment firstAlignment : candidates(first))
        {
            for(LiftedAlignment secondAlignment : candidates(second))
            {
                PlacementPairPriority priority = priority(firstAlignment, secondAlignment);
                PairCandidate pair = new PairCandidate(firstAlignment, secondAlignment);
                int comparison = bestPriority == null ? -1 : priority.compareTo(bestPriority);
                boolean better = comparison < 0;
                boolean tied = comparison == 0;
                if(better)
                {
                    bestPairs.clear();
                    bestPriority = priority;
                }
                if(better || tied)
                {
                    bestPairs.add(pair);
                }
            }
        }

        if(bestPairs.isEmpty())
        {
            return Optional.empty();
        }

        List<PairCandidate> contenders = topScoring(bestPairs);

        contenders.sort(Comparator.comparing(PairCandidate::canonicalKey));
        PairCandidate winner = contenders.get(Math.floorMod(seed, contenders.size()));
        String note = bestPairs.size() == 1 ? "mate" : (contenders.size() == 1 ? "score" : "random");
        return Optional.of(new Choice(winner.first(), winner.second(), note));
    }

    private static List<LiftedAlignment> candidates(final MatePlacements mate)
    {
        if(!mate.selectable())
        {
            return List.of(mate.selected());
        }

        List<LiftedAlignment> candidates = new ArrayList<>();
        Set<AlignmentKey> seen = new HashSet<>();
        for(LiftedAlignment alignment : mate.alignments())
        {
            if(!alignment.Dropped
                    && (alignment.GenomicScore == Integer.MIN_VALUE
                            || alignment.GenomicScore >= PRIMARY_AS_UNMAP_THRESHOLD)
                    && seen.add(alignment.key()))
            {
                candidates.add(alignment);
            }
        }
        if(!candidates.isEmpty())
        {
            return candidates;
        }

        // Keep the normal unmap path reachable when every scored placement is below the AS floor.
        for(LiftedAlignment alignment : mate.alignments())
        {
            if(!alignment.Dropped && seen.add(alignment.key()))
            {
                candidates.add(alignment);
            }
        }
        return candidates;
    }

    // pair priority ranks by SV category and separation, so where it cannot separate two pairs the better-scoring
    // placement is taken rather than an arbitrary one
    private static List<PairCandidate> topScoring(final List<PairCandidate> pairs)
    {
        long topScore = Long.MIN_VALUE;
        for(PairCandidate pair : pairs)
        {
            if(pair.first().GenomicScore == Integer.MIN_VALUE || pair.second().GenomicScore == Integer.MIN_VALUE)
            {
                return pairs;
            }
            topScore = Math.max(topScore, (long) pair.first().GenomicScore + pair.second().GenomicScore);
        }

        List<PairCandidate> topScored = new ArrayList<>();
        for(PairCandidate pair : pairs)
        {
            if((long) pair.first().GenomicScore + pair.second().GenomicScore == topScore)
            {
                topScored.add(pair);
            }
        }
        return topScored;
    }

    private static PlacementPairPriority priority(
            final LiftedAlignment first, final LiftedAlignment second)
    {
        int firstBreakend = first.ForwardStrand ? first.alignedEnd() : first.LiftedPos;
        int secondBreakend = second.ForwardStrand ? second.alignedEnd() : second.LiftedPos;
        Orientation firstOrientation = first.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
        Orientation secondOrientation = second.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;

        // an intron is not aligned sequence, so separation is measured between aligned blocks rather than across the
        // outer cigar span, which would penalise a spliced placement by the length of its introns
        int separation = first.alignedBlockDistance(second);

        return PlacementPairPriority.betweenMates(
                first.LiftedChromosome, firstBreakend, firstOrientation,
                second.LiftedChromosome, secondBreakend, secondOrientation, separation);
    }

    record MatePlacements(
            List<LiftedAlignment> alignments, LiftedAlignment selected, boolean selectable)
    {
        MatePlacements
        {
            alignments = List.copyOf(alignments);
        }
    }

    record Choice(LiftedAlignment first, LiftedAlignment second, String note) { }

    private record PairCandidate(LiftedAlignment first, LiftedAlignment second)
    {
        String canonicalKey()
        {
            String firstKey = first.key().toString();
            String secondKey = second.key().toString();
            return firstKey.compareTo(secondKey) <= 0
                    ? firstKey + '|' + secondKey
                    : secondKey + '|' + firstKey;
        }
    }
}
