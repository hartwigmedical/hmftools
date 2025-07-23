package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.BiPredicate;
import java.util.function.DoublePredicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

// TODO: unit test

// Miscellaneous utility functionality.
public class Utils
{
    // Min/max function with early stopping if an optimal value is found.
    public static <T> Optional<T> getBestScoringElement(Stream<T> elements, final ToDoubleFunction<T> scoreFunc,
            final DoublePredicate isOptimalFunc, boolean maximise)
    {
        BiPredicate<Double, Double> scoreCompareFunc = maximise ? Doubles::greaterThan : Doubles::lessThan;
        Optional<T> bestElement = Optional.empty();
        double bestScore = maximise ? NEGATIVE_INFINITY : POSITIVE_INFINITY;
        Iterator<T> iterator = elements.iterator();
        while(iterator.hasNext())
        {
            T element = iterator.next();
            double score = scoreFunc.applyAsDouble(element);
            if(scoreCompareFunc.test(score, bestScore))
            {
                bestElement = Optional.of(element);
                bestScore = score;
                if(isOptimalFunc.test(bestScore))
                {
                    break;
                }
            }
        }
        return bestElement;
    }

    // Compute regions within `targetRegion` which do not overlap `coveredRegions`.
    public static List<BaseRegion> computeUncoveredRegions(final BaseRegion targetRegion, Stream<BaseRegion> coveredRegions)
    {
        if(!targetRegion.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        coveredRegions = coveredRegions.peek(coveredRegion ->
        {
            if(!coveredRegion.hasValidPositions())
            {
                throw new IllegalArgumentException("Invalid region");
            }
        });

        // Ignore covered positions which don't overlap the target region, since they can never produce an uncovered region.
        coveredRegions = coveredRegions.filter(targetRegion::overlaps);

        // Sort by start position ascending, then end position ascending.
        coveredRegions = coveredRegions.sorted(Comparator.comparing(BaseRegion::start).thenComparing(BaseRegion::end));

        List<BaseRegion> uncoveredRegions = new ArrayList<>();

        Iterator<BaseRegion> iterator = coveredRegions.iterator();
        // Setting this to just before the target region simplifies handling of first and last uncovered regions.
        int prevCoveredPos = targetRegion.start() - 1;

        // Remaining covered regions.
        while(iterator.hasNext())
        {
            BaseRegion coveredRegion = iterator.next();
            // Possibilities:
            //   - Current region starts at same position as previous:
            //     - And ends >= previous end: nothing to do.
            //   - Current region starts after previous start:
            //     - And ends <= previous end + 1: nothing to do.
            //     - And ends > previous end + 1: uncovered region in between.
            if(coveredRegion.start() > prevCoveredPos + 1)
            {
                int uncoveredStart = prevCoveredPos + 1;
                int uncoveredEnd = min(coveredRegion.start() - 1, targetRegion.end());
                uncoveredRegions.add(new BaseRegion(uncoveredStart, uncoveredEnd));
            }
            prevCoveredPos = max(prevCoveredPos, coveredRegion.end());
        }

        if(prevCoveredPos < targetRegion.end())
        {
            // Covered regions end before the end of the target region, so there is an uncovered region afterward.
            uncoveredRegions.add(new BaseRegion(prevCoveredPos + 1, targetRegion.end()));
        }

        return uncoveredRegions;
    }

    // Checks if `targetRegion` is fully covered by `coveredRegions`.
    public static boolean isCoveredBy(final ChrBaseRegion targetRegion, Stream<ChrBaseRegion> coveredRegions)
    {
        // Similar to computeUncoveredRegions() but without tracking the uncovered regions.
        coveredRegions = coveredRegions.filter(targetRegion::overlaps);
        coveredRegions = coveredRegions.sorted(Comparator.comparing(ChrBaseRegion::start).thenComparing(ChrBaseRegion::end));
        Iterator<ChrBaseRegion> iterator = coveredRegions.iterator();
        int prevCovered = targetRegion.start() - 1;
        if(!iterator.hasNext())
        {
            // No overlapping regions, therefore can't be covered.
            return false;
        }
        while(iterator.hasNext())
        {
            ChrBaseRegion coveredRegion = iterator.next();
            if(coveredRegion.start() > prevCovered + 1)
            {
                // Since the regions are sorted by start, a gap must mean an uncovered region.
                return false;
            }
            prevCovered = max(prevCovered, coveredRegion.end());
        }
        // If we got here then at least 1 region overlaps and there are no gaps in the coverage.
        return true;
    }

    public static double regionCentreFloat(final BaseRegion region)
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return (region.start() + region.end()) / 2.0;
    }

    public static int regionCentre(final BaseRegion region)
    {
        return (int) regionCentreFloat(region);
    }

    public static BaseRegion regionStartingAt(int startPosition, int length)
    {
        BaseRegion region = new BaseRegion(startPosition, startPosition + length - 1);
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return region;
    }

    public static int regionCentreStartOffset(int length)
    {
        return -(length / 2) + (1 - length % 2);
    }

    public static BaseRegion regionCenteredAt(int centrePosition, int length)
    {
        int start = centrePosition + regionCentreStartOffset(length);
        int end = start + length - 1;
        BaseRegion region = new BaseRegion(start, end);
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return region;
    }

    public static BaseRegion regionEndingAt(int endPosition, int length)
    {
        BaseRegion region = new BaseRegion(endPosition - length + 1, endPosition);
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return region;
    }

    public static BaseRegion regionIntersection(final BaseRegion region1, final BaseRegion region2)
    {
        if(!region1.hasValidPositions() || !region2.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return new BaseRegion(max(region1.start(), region2.start()), min(region1.end(), region2.end()));
    }

    public static boolean regionOverlapsOrAdjacent(final BaseRegion region1, final BaseRegion region2)
    {
        return region1.overlaps(region2) || region1.end() + 1 == region2.start() || region2.end() + 1 == region1.start();
    }

    // Generates the sequence: 0, 1, -1, 2, -2, 3, -3, ...
    // With the constraint that no value will be outside the range [minOffset, maxOffset].
    public static IntStream outwardMovingOffsets(int minOffset, int maxOffset)
    {
        if(minOffset > maxOffset)
        {
            // Probably a bug in the calling code.
            throw new IllegalArgumentException("minOffset and maxOffset forbid all possible offsets");
        }
        return IntStream.iterate(0, absOffset -> -absOffset >= minOffset || absOffset <= maxOffset, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset))
                .filter(offset -> offset >= minOffset && offset <= maxOffset);
    }
}
