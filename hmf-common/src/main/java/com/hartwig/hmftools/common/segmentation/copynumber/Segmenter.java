package com.hartwig.hmftools.common.segmentation.copynumber;

import static java.util.Collections.singletonList;

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

    public Segmenter(double[] y)
    {
        this(y, 50.0, false);
    }

    public Segmenter(double[] y, double gamma, boolean normalise)
    {
        this(y, new GammaPenaltyCalculator(gamma, normalise));
    }

    public Segmenter(double[] y, PenaltyCalculator penaltyCalculator)
    {
        this.y = y;
        segmentPenalty = penaltyCalculator.getPenalty(y);

        // Here we calculate leastCost(s,e) for each e in [0, ..., y.size - 1]
        // and s in [0, ..., e]. For any such pair, the least cost is:
        // leastCost(s,e) = leastCostOfSegmentingTo(s-1) + segmentPenalty + intervalCost(s,e)
        // The least cost segmentation corresponds to the leastCostOfSegmentingTo(y.size - 1)
        double[] leastCostEndingJustBefore = new double[y.length + 1];
        leastCostEndingJustBefore[0] = 0.0;
        // We record the segment endpoints so that we can recover the least cost segmentation itself.
        int[] leastCostSegmentEndpoints = new int[y.length]; // there's only one possible segment of length 1

        double[] cumulativeSums = precomputeCumulativeSums();
        double[] cumulativeSquaredSums = precomputeCumulativeSquaredSums();

        for(int end = 0; end < y.length; end++)
        {
            double minCost = Double.MAX_VALUE;
            int endOfPreviousSegmentForLeastCost = 0;

            for(int start = 0; start <= end; start++)
            {
                // Calculate segment cost using precomputed squared sums
                double sumSquared = start == 0
                        ? cumulativeSums[end] * cumulativeSums[end]
                        : (cumulativeSums[end] - cumulativeSums[start - 1]) * (cumulativeSums[end] - cumulativeSums[start - 1]);
                double sumOfSquares = start == 0 ?
                        cumulativeSquaredSums[end] :
                        cumulativeSquaredSums[end] - cumulativeSquaredSums[start - 1];
                int segmentLength = end - start + 1;
                double segmentCost = sumOfSquares - (sumSquared / segmentLength);

                double cost = leastCostEndingJustBefore[start] + segmentPenalty + segmentCost;
                if(cost < minCost)
                {
                    minCost = cost;
                    endOfPreviousSegmentForLeastCost = start - 1;
                }
            }
            leastCostEndingJustBefore[end + 1] = minCost;
            leastCostSegmentEndpoints[end] = endOfPreviousSegmentForLeastCost;
        }
        List<Integer> segmentEndpoints = new ArrayList<>();
        int lastSegmentEndpoint = y.length - 1;
        while(lastSegmentEndpoint >= 0)
        {
            segmentEndpoints.add(lastSegmentEndpoint);
            lastSegmentEndpoint = leastCostSegmentEndpoints[lastSegmentEndpoint];
        }
        Collections.reverse(segmentEndpoints);
        leastCostSegmentation = segmentEndpoints.isEmpty() ? new Segmentation(singletonList(y)) : segmentBy(segmentEndpoints);
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

        int prev = list.get(0);
        for(int i = 1; i < list.size(); i++)
        {
            int current = list.get(i);
            if(current <= prev)
            {
                return false;
            }
            prev = current;
        }
        return true;
    }

    private double[] precomputeCumulativeSums()
    {
        double[] cumulativeSums = new double[y.length];
        if(y.length > 0)
        {
            cumulativeSums[0] = y[0];
            for(int i = 1; i < y.length; i++)
            {
                cumulativeSums[i] = cumulativeSums[i - 1] + y[i];
            }
        }
        return cumulativeSums;
    }

    private double[] precomputeCumulativeSquaredSums()
    {
        double[] cumulativeSquaredSums = new double[y.length];
        if(y.length > 0)
        {
            cumulativeSquaredSums[0] = y[0] * y[0];
            for(int i = 1; i < y.length; i++)
            {
                cumulativeSquaredSums[i] = cumulativeSquaredSums[i - 1] + (y[i] * y[i]);
            }
        }
        return cumulativeSquaredSums;
    }
}
