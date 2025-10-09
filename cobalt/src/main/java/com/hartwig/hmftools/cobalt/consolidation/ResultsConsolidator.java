package com.hartwig.hmftools.cobalt.consolidation;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Comparators;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.BamRatios;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.apache.commons.lang3.Validate;

public interface ResultsConsolidator
{
    static int calcConsolidationCount(final double medianReadDepth)
    {
        // consolidation starts when mean read depth <= 8
        double c = 80.0 / medianReadDepth;

        if(c < 10.0)
        {
            return 1;
        }

        // max 1000
        if(c >= 1_000)
        {
            return 1_000;
        }

        // round to one significant digit
        double roundBy = Math.pow(10, Math.floor(Math.log10(c)));
        return (int) (Math.round(c / roundBy) * roundBy);
    }

    // given the list of non masked windows, get the list of consolidated buckets
    static List<LowCovBucket> consolidateIntoBuckets(List<Integer> windowPositions, int consolidationCount)
    {
        // make sure position is sorted
        Validate.isTrue(Comparators.isInStrictOrder(windowPositions, Comparator.naturalOrder()));

        List<LowCovBucket> buckets = new ArrayList<>();

        if(windowPositions.isEmpty())
        {
            return buckets;
        }

        int windowCount = 0;
        int bucketStart = windowPositions.get(0);

        for(int i = 0; i < windowPositions.size(); ++i)
        {
            int position = windowPositions.get(i);

            if((position - bucketStart) >= CobaltConstants.MAX_SPARSE_CONSOLIDATE_DISTANCE)
            {
                // do not let the bucket to consolidate over 3M bases. This is done to avoid consolidating
                // over centromere
                // use the last bucket
                if(i > 0)
                {
                    int lastPosition = windowPositions.get(i - 1);
                    int bucketEnd = lastPosition + CobaltConstants.WINDOW_SIZE;
                    int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
                    buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));
                }

                // reset bucket start position, we do not want to put the start in the middle since
                // it would sit inside the centromere
                bucketStart = position;

                // also reset window count
                windowCount = 0;
            }

            if(windowCount == consolidationCount)
            {
                // we want to put the bucket boundary in the middle of the two windows
                int lastPosition = windowPositions.get(i - 1);
                int bucketEnd = roundDownToWindowBoundary((lastPosition + position) * 0.5);

                // bucket position is at the middle
                int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);

                buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));

                // next bucket starts right after this
                bucketStart = bucketEnd + CobaltConstants.WINDOW_SIZE;

                // also reset window count
                windowCount = 0;
            }

            windowCount++;
        }

        // add a final window
        if(windowCount > 0)
        {
            int bucketEnd = windowPositions.get(windowPositions.size() - 1) + CobaltConstants.WINDOW_SIZE;

            // bucket position is at the middle
            int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
            buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));
        }

        return buckets;
    }

    default ListMultimap<Chromosome, BamRatio> consolidate(ListMultimap<Chromosome, BamRatio> ratios)
    {
        return ratios;
    }

    private static int roundDownToWindowBoundary(double p)
    {
        return (int) (Math.floor(p / CobaltConstants.WINDOW_SIZE) * CobaltConstants.WINDOW_SIZE) + 1;
    }
}
