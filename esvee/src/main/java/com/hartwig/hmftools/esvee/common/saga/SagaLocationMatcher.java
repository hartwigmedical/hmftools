package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.utils.IntPair;

import org.jetbrains.annotations.Nullable;

// Matches junctions to SAGA variants by genomic location.
// For a junction, find the nearest SAGA variant breakend within the configured distance tolerance.
public class SagaLocationMatcher
{
    private final Config mConfig;
    // Map from contig to list of breakends sorted by position ascending.
    private final Map<String, List<SagaIndexedBreakend>> mSearchableBreakends;

    // Do not instantiate directly - use SagaMatcherFactory.
    // searchableBreakends MUST be sorted ascending by breakend position.
    public SagaLocationMatcher(final Config config, final Map<String, List<SagaIndexedBreakend>> searchableBreakends)
    {
        mConfig = config;
        mSearchableBreakends = searchableBreakends;
    }

    public record Config(
            int locationDistanceMax
    )
    {
        public static final Config DEFAULT = new Config(SAGA_LOCATION_MATCH_DISTANCE);
    }

    @Nullable
    public SagaLocationMatch match(final String chromosome, int position, final Orientation orientation)
    {
        Map<String, List<SagaIndexedBreakend>> breakends = mSearchableBreakends;
        List<SagaIndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchOnChromosome(chrBreakends, position, orientation);
    }

    @Nullable
    private SagaLocationMatch matchOnChromosome(@Nullable final List<SagaIndexedBreakend> chrBreakends, int position,
            final Orientation orientation)
    {
        if(chrBreakends == null || chrBreakends.isEmpty())
        {
            return null;
        }

        IntPair range = SagaIndexedBreakend.binarySearch(chrBreakends, position);
        if(range.left == range.right)
        {
            // Position isn't found exactly, so need to look at adjacent indices.
            int startIndex = max(range.left - 1, 0);
            int endIndex = min(range.left + 1, chrBreakends.size());
            return chooseFromInexactMatchCandidates(chrBreakends.subList(startIndex, endIndex), position, orientation);
        }
        else
        {
            // Position found exactly. There may be multiple breakends with the same position.
            return chooseFromExactMatchCandidates(chrBreakends.subList(range.left, range.right), orientation);
        }
    }

    @Nullable
    private SagaLocationMatch chooseFromInexactMatchCandidates(final List<SagaIndexedBreakend> candidateBreakends, int position,
            final Orientation orientation)
    {
        // Consider the nearby breakends, compute their distances, and filter out the breakends too far away.
        List<SagaLocationMatch> allCandidates = candidateBreakends.stream()
                .map(breakend -> new SagaLocationMatch(breakend.variant(), breakend.breakend(), abs(breakend.position() - position)))
                .filter(match -> match.distance() <= mConfig.locationDistanceMax)
                .toList();
        // Select only the breakend(s) with the lowest distance. This may be 1 breakend or multiple.
        int minDistance = allCandidates.stream().mapToInt(SagaLocationMatch::distance).min().orElse(mConfig.locationDistanceMax + 1);
        if(minDistance > mConfig.locationDistanceMax)
        {
            return null;
        }
        Stream<SagaLocationMatch> equidistantCandidates = allCandidates.stream().filter(match -> match.distance() == minDistance);
        // Select the match based on tie-break conditions.
        return chooseFromEquidistantCandidates(equidistantCandidates, orientation);
    }

    private static SagaLocationMatch chooseFromExactMatchCandidates(final List<SagaIndexedBreakend> candidateBreakends,
            final Orientation orientation)
    {
        Stream<SagaLocationMatch> candidates = candidateBreakends.stream()
                .map(breakend -> new SagaLocationMatch(breakend.variant(), breakend.breakend(), 0));
        return chooseFromEquidistantCandidates(candidates, orientation);
    }

    private static SagaLocationMatch chooseFromEquidistantCandidates(final Stream<SagaLocationMatch> candidates,
            final Orientation orientation)
    {
        // If there are multiple candidates with the same lowest distance, choose the breakend with the same orientation.
        // Then finally order by variant ID to deterministically resolve further ties.
        Comparator<SagaLocationMatch> comparator = Comparator
                .comparing((SagaLocationMatch match) -> match.breakend().orientation() != orientation)
                .thenComparing(SagaLocationMatch::variantId);
        return candidates.min(comparator).orElseThrow();
    }
}
