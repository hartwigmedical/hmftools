package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;

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

        List<BaseRegion> uncoveredRegions = new ArrayList<>();

        Iterator<BaseRegion> iterator = coveredRegions.iterator();
        BaseRegion prevCoveredRegion = null;
        while (iterator.hasNext()) {
            BaseRegion coveredRegion = iterator.next();
            if (prevCoveredRegion == null) {
                // First covered region, check against target region start.
                if (coveredRegion.start() > targetRegion.start()) {
                    uncoveredRegions.add(new BaseRegion(targetRegion.start(), coveredRegion.start() - 1));
                }
            }
            else {
                // Possibilities:
                //   - Current region starts at same position as previous:
                //     - And ends >= previous end: nothing to do.
                //   - Current region starts after previous start:
                //     - And ends <= previous end + 1: nothing to do.
                //     - And ends > previous end + 1: uncovered region in between.
                if (coveredRegion.start() > prevCoveredRegion.end() + 1) {
                    uncoveredRegions.add(new BaseRegion(prevCoveredRegion.end() + 1, coveredRegion.start() - 1));
                }
            }
            prevCoveredRegion = coveredRegion;
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
