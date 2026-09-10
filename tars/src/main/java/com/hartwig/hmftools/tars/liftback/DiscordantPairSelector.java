package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.tars.liftback.features.LocalSvPriority;

// Selects the two emitted placements as one fragment when either mate has MAPQ 0.
final class DiscordantPairSelector
{
    private DiscordantPairSelector() { }

    static Optional<Choice> select(final MatePlacements first, final MatePlacements second, final int seed)
    {
        List<PairCandidate> bestPairs = new ArrayList<>();
        LocalSvPriority.Rank bestPriority = null;

        for(LiftedAlignment firstAlignment : candidates(first))
        {
            for(LiftedAlignment secondAlignment : candidates(second))
            {
                LocalSvPriority.Rank priority = priority(firstAlignment, secondAlignment);
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

        bestPairs.sort(Comparator.comparing(PairCandidate::canonicalKey));
        PairCandidate winner = bestPairs.get(Math.floorMod(seed, bestPairs.size()));
        return Optional.of(new Choice(
                winner.first(), winner.second(), bestPairs.size() == 1 ? "mate" : "random"));
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

    private static LocalSvPriority.Rank priority(
            final LiftedAlignment first, final LiftedAlignment second)
    {
        int firstBreakend = first.ForwardStrand ? first.alignedEnd() : first.LiftedPos;
        int secondBreakend = second.ForwardStrand ? second.alignedEnd() : second.LiftedPos;
        Orientation firstOrientation = first.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
        Orientation secondOrientation = second.ForwardStrand ? Orientation.FORWARD : Orientation.REVERSE;
        return LocalSvPriority.between(
                first.LiftedChromosome, firstBreakend, firstOrientation,
                second.LiftedChromosome, secondBreakend, secondOrientation);
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
