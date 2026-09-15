package com.hartwig.hmftools.purple.reseg;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalInt;
import java.util.Set;
import java.util.function.IntPredicate;
import java.util.function.ToDoubleFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries.Bucket;

public final class PeakTroughFinder
{
    private PeakTroughFinder() {}

    public static List<Integer> findLocalPeakOrTroughIndices(final List<Bucket> series, boolean findPeaks)
    {
        // find peak or trough indices
        List<Integer> indices = Lists.newArrayList();

        for(int i = 1; i < series.size() - 1; i++)
        {
            double prev = series.get(i - 1).value();
            double curr = series.get(i).value();
            double next = series.get(i + 1).value();

            if(findPeaks)
            {
                if(curr > prev && curr > next)
                    indices.add(i);
            }
            else
            {
                if(curr < prev && curr < next)
                    indices.add(i);
            }
        }

        return indices;
    }

    public static OptionalInt findNearestBelowRequired(
            final List<Integer> sortedCandidates, int centerIndex, final IntPredicate satisfies)
    {
        // searches sortedCandidates (ascending) for the nearest index below centerIndex satisfying the predicate,
        // trying progressively further candidates only if nearer ones fail
        for(int i = sortedCandidates.size() - 1; i >= 0; i--)
        {
            int candidate = sortedCandidates.get(i);

            if(candidate >= centerIndex)
                continue;

            if(satisfies.test(candidate))
                return OptionalInt.of(candidate);
        }

        return OptionalInt.empty();
    }

    public static OptionalInt findNearestAboveRequired(
            final List<Integer> sortedCandidates, int centerIndex, final IntPredicate satisfies)
    {
        for(int candidate : sortedCandidates)
        {
            if(candidate <= centerIndex)
                continue;

            if(satisfies.test(candidate))
                return OptionalInt.of(candidate);
        }

        return OptionalInt.empty();
    }

    public static List<PeakTroughData> consolidateResults(
            final List<PeakTroughData> sortedByLevel, final KeepBothTest<PeakTroughData> keepBothTest,
            final ToDoubleFunction<PeakTroughData> valueFn, boolean keepHigherValue)
    {
        // Pairwise left-to-right consolidation, mirroring the reference's "collect all drops from this pass over
        // the pre-pass list, then apply them all at once, repeat until a full pass drops nothing" semantics -
        // NOT immediate in-loop removal, which would let earlier drops influence later decisions within the same pass.
        List<PeakTroughData> current = new ArrayList<>(sortedByLevel);
        boolean anyDropped = true;

        while(anyDropped)
        {
            anyDropped = false;
            Set<Integer> dropPositions = new HashSet<>();

            for(int i = 0; i < current.size() - 1; i++)
            {
                PeakTroughData left = current.get(i);
                PeakTroughData right = current.get(i + 1);

                if(keepBothTest.keepBoth(left, right))
                    continue;

                anyDropped = true;

                double leftValue = valueFn.applyAsDouble(left);
                double rightValue = valueFn.applyAsDouble(right);

                boolean dropLeft = keepHigherValue ? (rightValue > leftValue) : (rightValue < leftValue);

                dropPositions.add(dropLeft ? i : i + 1);
            }

            if(anyDropped)
            {
                List<PeakTroughData> next = new ArrayList<>();

                for(int i = 0; i < current.size(); i++)
                {
                    if(!dropPositions.contains(i))
                        next.add(current.get(i));
                }

                current = next;
            }
        }

        return current;
    }

    public interface KeepBothTest<T>
    {
        boolean keepBoth(T left, T right);
    }
}
