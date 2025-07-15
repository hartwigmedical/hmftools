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
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

public class Utils
{
    // Compute regions within `targetRegion` which do not overlap `coveredRegions`.
    public static List<BaseRegion> computeUncoveredRegions(final BaseRegion targetRegion, Stream<BaseRegion> coveredRegions)
    {
        // Sort by start position ascending, then end position ascending.
        coveredRegions = coveredRegions.sorted(Comparator.comparing(BaseRegion::start).thenComparing(BaseRegion::end));

        // Ignore covered positions which don't overlap the target region, since they can never produce an uncovered region.
        coveredRegions = coveredRegions.filter(targetRegion::overlaps);

        List<BaseRegion> uncoveredRegions = new ArrayList<>();

        Iterator<BaseRegion> iterator = coveredRegions.iterator();
        // Setting this to just before the target region makes the code simpler for the case of no covered regions.
        int prevCoveredPos = targetRegion.start() - 1;

        // Special handling for first covered region which is checked against the target region start.
        if(iterator.hasNext())
        {
            BaseRegion coveredRegion = iterator.next();
            if(coveredRegion.start() > targetRegion.start())
            {
                int uncoveredStart = targetRegion.start();
                int uncoveredEnd = min(targetRegion.end(), coveredRegion.start() - 1);
                uncoveredRegions.add(new BaseRegion(uncoveredStart, uncoveredEnd));
            }
            prevCoveredPos = coveredRegion.end();
        }

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
                uncoveredRegions.add(new BaseRegion(prevCoveredPos + 1, coveredRegion.start() - 1));
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
}
