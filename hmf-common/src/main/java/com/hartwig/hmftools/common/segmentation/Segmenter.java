package com.hartwig.hmftools.common.segmentation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Segmenter
{
    private final double[] y;
    private final Segmentation leastCostSegmentation;
    private final double segmentPenalty;

    public Segmenter(double[] y, double gamma, boolean normalise)
    {
        this.y = y;

        // Here we calculate leastCost(s,e) for each e in [0,..,y.size - 1]
        // and s in [0,..,e]. For any such pair, the least cost is:
        // leastCost(s,e) = leastCostOfSegmentingTo(s-1) + segmentPenalty + intervalCost(s,e)
        // The least cost segmentation corresponds to the leastCostOfSegmentingTo(y.size - 1)
        Map<Integer, Double> leastCostEndingAt = new HashMap<>();
        leastCostEndingAt.put(-1, 0.0);
        // We record the segment endpoints so that we can recover the least cost segmentation itself.
        List<Integer> leastCostSegmentEndpoints = new ArrayList<>(); // there's only one possible segment of length 1
        segmentPenalty = new Gamma(y, gamma, normalise).getSegmentPenalty();
        System.out.println("Segment penalty: " + segmentPenalty);
        List<Double> partialSumsFromIndexToCurrentEnd = new ArrayList<>();
        for(int end = 0; end < y.length; end++)
        {
            for(int i = 0; i < partialSumsFromIndexToCurrentEnd.size(); i++)
            {
                partialSumsFromIndexToCurrentEnd.set(i, partialSumsFromIndexToCurrentEnd.get(i) + y[end]);
            }
            partialSumsFromIndexToCurrentEnd.add(y[end]);
            double minCost = Double.MAX_VALUE;
            int endOfPreviousSegmentForLeastCost = 0;
            for(int start = 0; start <= end; start++)
            {
                double partialSum = partialSumsFromIndexToCurrentEnd.get(start);
                double segmentCost = 1 * partialSum * partialSum / ((start - end - 1.0f));
                double cost = leastCostEndingAt.get(start - 1) + segmentPenalty + segmentCost;
                if(cost < minCost)
                {
                    minCost = cost;
                    endOfPreviousSegmentForLeastCost = start - 1;
                }
            }
            leastCostEndingAt.put(end, minCost);
            leastCostSegmentEndpoints.add(endOfPreviousSegmentForLeastCost);
        }
        List<Integer> segmentEndpoints = new ArrayList<>();
        int lastSegmentEndpoint = y.length - 1;
        while(lastSegmentEndpoint >= 0)
        {
            segmentEndpoints.add(lastSegmentEndpoint);
            lastSegmentEndpoint = leastCostSegmentEndpoints.get(lastSegmentEndpoint);
        }
        Collections.reverse(segmentEndpoints);
        leastCostSegmentation = segmentEndpoints.isEmpty() ?
                new Segmentation(Collections.singletonList(y)) :
                segmentBy(segmentEndpoints);
    }

    public Segmenter(double[] y)
    {
        this(y, 50.0, false);
    }

    public PCF pcf()
    {
        return leastCostSegmentation.pcf();
    }

    public Segmentation cheapestSegmentationByExhaustiveSearch()
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

    public PCF pcfByExhaustiveSearch()
    {
        return cheapestSegmentationByExhaustiveSearch().pcf();
    }

    public Set<Segmentation> allPossibleSegmentations()
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

    public Segmentation segmentBy(List<Integer> segmentEndpointIndices)
    {
        if(!isIncreasing(segmentEndpointIndices))
        {
            throw new IllegalArgumentException("Segment endpoint indices must be increasing");
        }
        if(segmentEndpointIndices.get(segmentEndpointIndices.size() - 1) != y.length - 1)
        {
            throw new IllegalArgumentException("Last segment endpoint index must be the last index of the input array");
        }

        List<double[]> segments = new ArrayList<>();
        int start = 0;
        for(int it : segmentEndpointIndices)
        {
            int segmentSize = it - start + 1;
            double[] da = new double[segmentSize];
            if(it + 1 - start >= 0)
            {
                System.arraycopy(y, start, da, 0, it + 1 - start);
            }
            segments.add(da);
            start = it + 1;
        }
        return new Segmentation(segments);
    }

    private static boolean isIncreasing(List<Integer> list)
    {
        if(list.size() < 2)
        {
            return true;
        }
        for(int i = list.size() - 2; i >= 1; i--)
        {
            if(list.get(i) <= list.get(i - 1))
            {
                return false;
            }
        }
        return true;
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
        List<double[]> headSplits = splitByIndexes(headTail.getFirst(), remainingIndexes);
        List<double[]> result = new ArrayList<>(headSplits);
        result.add(headTail.getSecond());
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
        return new Pair<>(left, right);
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

    private static class Pair<F, S>
    {
        private final F first;
        private final S second;

        public Pair(F first, S second)
        {
            this.first = first;
            this.second = second;
        }

        public F getFirst()
        {
            return first;
        }

        public S getSecond()
        {
            return second;
        }
    }
}