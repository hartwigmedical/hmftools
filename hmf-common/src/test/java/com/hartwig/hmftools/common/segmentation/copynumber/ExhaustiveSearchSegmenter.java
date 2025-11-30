package com.hartwig.hmftools.common.segmentation.copynumber;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;

public class ExhaustiveSearchSegmenter
{
    private final double[] y;
    public final double segmentPenalty;

    public ExhaustiveSearchSegmenter(double[] y, double gamma, boolean normalise)
    {
        this.y = y;
        segmentPenalty = new Gamma(y, gamma, normalise).getSegmentPenalty();
    }

    public ExhaustiveSearchSegmenter(double[] y)
    {
        this(y, 50.0, false);
    }

    Segmentation cheapestSegmentationByExhaustiveSearch()
    {
        Set<Segmentation> allSegmentations = allPossibleSegmentations();
        Segmentation cheapest = null;
        double minCost = Double.MAX_VALUE;
        for(Segmentation segmentation : allSegmentations)
        {
            double cost = segmentation.cost(segmentPenalty);
            if(cost < minCost)
            {
                minCost = cost;
                cheapest = segmentation;
            }
        }
        return cheapest;
    }

    Set<Segmentation> allPossibleSegmentations()
    {
        if(y.length == 1)
        {
            return Collections.singleton(new Segmentation(Collections.singletonList(y)));
        }
        Set<Segmentation> result = new HashSet<>();
        Set<Integer> indexes = new HashSet<>();
        for(int i = 1; i < y.length; i++)
        {
            indexes.add(i);
        }
        Set<Set<Integer>> power = powerSet(indexes);
        for(Set<Integer> indexSet : power)
        {
            result.add(new Segmentation(splitByIndexes(y, indexSet)));
        }
        return result;
    }

    private static List<double[]> splitByIndexes(double[] array, Set<Integer> indexes)
    {
        for(int index : indexes)
        {
            if(index <= 0 || index >= array.length)
            {
                throw new IllegalArgumentException("Index must be > 0 and < array.length");
            }
        }
        if(indexes.isEmpty())
        {
            return Collections.singletonList(array);
        }

        int greatestIndex = Collections.max(indexes);
        Pair<double[], double[]> headTail = splitAt(array, greatestIndex);
        Set<Integer> remainingIndexes = new HashSet<>(indexes);
        remainingIndexes.remove(greatestIndex);
        List<double[]> headSplits = splitByIndexes(headTail.getLeft(), remainingIndexes);
        List<double[]> result = new ArrayList<>(headSplits);
        result.add(headTail.getRight());
        return result;
    }

    private static Pair<double[], double[]> splitAt(double[] array, int index)
    {
        if(index <= 0 || index >= array.length)
        {
            throw new IllegalArgumentException("Index must be > 0 and < array.length");
        }
        double[] left = new double[index];
        double[] right = new double[array.length - index];
        for(int i = 0; i < array.length; i++)
        {
            if(i < index)
            {
                left[i] = array[i];
            }
            else
            {
                right[i - index] = array[i];
            }
        }
        return Pair.of(left, right);
    }

    private static <T> Set<Set<T>> powerSet(Set<T> set)
    {
        Set<Set<T>> result = new HashSet<>();
        result.add(new HashSet<>());
        for(T element : set)
        {
            Set<T> withoutElement = new HashSet<>(set);
            withoutElement.remove(element);
            Set<Set<T>> power = powerSet(withoutElement);
            for(Set<T> subset : power)
            {
                result.add(subset);
                Set<T> withElement = new HashSet<>(subset);
                withElement.add(element);
                result.add(withElement);
            }
        }
        return result;
    }

}
