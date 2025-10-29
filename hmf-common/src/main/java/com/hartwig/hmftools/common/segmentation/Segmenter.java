package com.hartwig.hmftools.common.segmentation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Given a list of doubles and a segmentation penalty this class finds the least cost
 * segmentation of the doubles using a dynamic programming algorithm described in
 * "Copynumber: Eﬃcient algorithms for single- and multi-track copy number segmentation"
 * Nilsen et al. BMC Genomics 2012,13:591 <a href="http://www.biomedcentral.com/1471-2164/13/591">...</a>
 */
public class Segmenter
{
    private final double[] y;
    final Segmentation leastCostSegmentation;
    public final double segmentPenalty;

    public Segmenter(double[] y, double gamma, boolean normalise)
    {
        this.y = y;

        // Here we calculate leastCost(s,e) for each e in [0, ..., y.size - 1]
        // and s in [0, ..., e]. For any such pair, the least cost is:
        // leastCost(s,e) = leastCostOfSegmentingTo(s-1) + segmentPenalty + intervalCost(s,e)
        // The least cost segmentation corresponds to the leastCostOfSegmentingTo(y.size - 1)
        List<Double> leastCostEndingJustBefore = new ArrayList<>(y.length);
        leastCostEndingJustBefore.add(0, 0.0);
        // We record the segment endpoints so that we can recover the least cost segmentation itself.
        List<Integer> leastCostSegmentEndpoints = new ArrayList<>(); // there's only one possible segment of length 1
        segmentPenalty = new Gamma(y, gamma, normalise).getSegmentPenalty();

        // Precompute cumulative sums for efficient partial sum calculation.
        double[] cumulativeSums = new double[y.length];
        if (y.length > 0) {
            cumulativeSums[0] = y[0];
            for (int i = 1; i < y.length; i++) {
                cumulativeSums[i] = cumulativeSums[i - 1] + y[i];
            }
        }

        for(int end = 0; end < y.length; end++)
        {
            double minCost = Double.MAX_VALUE;
            int endOfPreviousSegmentForLeastCost = 0;
            for(int start = 0; start <= end; start++)
            {
                // Calculate partial sum efficiently using precomputed cumulative sums
                double partialSum = start == 0 ? cumulativeSums[end] : cumulativeSums[end] - cumulativeSums[start - 1];
                double segmentCost = 1 * partialSum * partialSum / ((start - end - 1.0f));
                double cost = leastCostEndingJustBefore.get(start) + segmentPenalty + segmentCost;
                if(cost < minCost)
                {
                    minCost = cost;
                    endOfPreviousSegmentForLeastCost = start - 1;
                }
            }
            leastCostEndingJustBefore.add(end + 1, minCost);
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

    public PiecewiseConstantFit pcf()
    {
        return leastCostSegmentation.pcf();
    }

    Segmentation segmentBy(List<Integer> segmentEndpointIndices)
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
}
