package com.hartwig.hmftools.panelbuilder;

import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.max;
import static java.util.Objects.requireNonNull;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.BiPredicate;
import java.util.function.DoublePredicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.Doubles;

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

    // Generates the sequence: 0, 1, -1, 2, -2, 3, -3, ...
    // With the constraint that no value will be outside the range [minOffset, maxOffset].
    public static IntStream outwardMovingOffsets(int minOffset, int maxOffset)
    {
        return IntStream.iterate(0, absOffset -> -absOffset >= minOffset || absOffset <= maxOffset, absOffset -> absOffset + 1)
                .flatMap(absOffset -> absOffset == 0 ? IntStream.of(absOffset) : IntStream.of(absOffset, -absOffset))
                .filter(offset -> offset >= minOffset && offset <= maxOffset);
    }

    public static <T> List<T> findDuplicates(final List<T> items, BiPredicate<T, T> isDuplicate)
    {
        List<T> duplicates = new ArrayList<>();
        for(int i = 0; i < items.size(); i++)
        {
            T item1 = items.get(i);
            for(int j = 0; j < items.size(); j++)
            {
                if(i != j)
                {
                    T item2 = items.get(j);
                    if(isDuplicate.test(item1, item2))
                    {
                        if(!duplicates.contains(item1))
                        {
                            duplicates.add(item1);
                        }
                        if(!duplicates.contains(item2))
                        {
                            duplicates.add(item2);
                        }
                    }
                }
            }
        }
        return duplicates;
    }

    public static <T> List<T> findDuplicates(final List<T> items)
    {
        return findDuplicates(items, Object::equals);
    }

    public static double estimatePanelOnTargetRate(final List<Probe> probes)
    {
        double offTargetSum = probes.stream().mapToDouble(p -> qualityScoreToOffTargets(requireNonNull(p.qualityScore()))).sum();
        double offTargetRate = offTargetSum / probes.size();
        return 1 - offTargetRate;
    }

    private static double qualityScoreToOffTargets(double qualityScore)
    {
        if (qualityScore < 0)
        {
            throw new IllegalArgumentException("Quality score must be non-negative");
        }
        // The probe quality profile resolution is 0.01, so let's estimate that any probe which scores near 0 is about 0.01.
        // This avoids division by 0 issues. Because technically 0 QS means the probe has infinite off-targets.
        if (qualityScore < 1e-3)
        {
            qualityScore = 0.01;
        }
        return (1 / qualityScore) - 1;
    }
}
