package com.hartwig.hmftools.common.sequence;

import org.apache.commons.lang3.Range;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LongestSegment
{
    // longest segment algorithm where the average value of the segment
    // is >= averageL
    // this can be reduced to the longest segment where the sum is <= 0
    //
    // @return Range of the segment if found, null otherwise
    //
    @Nullable
    public static Range<Integer> longestSegmentAverage(@NotNull double[] a, double averageL)
    {
        double[] aMod = new double[a.length];
        for (int i = 0; i < aMod.length; ++i)
        {
            aMod[i] = a[i] - averageL;
        }
        return longestSegmentSum(aMod, 0.0);
    }

    // longest segment algorithm
    //
    // K.-Y. Chen and K.-M. Chao. Optimal algorithms for locating the
    // longest or shortest segments satisfying a sum or average constraint.
    // Inform. Process. Lett., 96:197–201, 2005.
    //
    // given a sequence of values, we want to find the longest segment
    // within this sequence where the sum is >= l
    //
    // note that I have shifted the indices such that they are 0 based instead
    //
    // @return Range of the segment if found, null otherwise
    //
    @Nullable
    public static Range<Integer> longestSegmentSum(@NotNull double[] a, double sumL)
    {
        double[] c = new double[a.length + 1];
        int[] m = new int[a.length + 1];
        int x = 0;
        int y = 0;

        for (int i = 1; i <= a.length; ++i)
        {
            // update the cumulative sum
            c[i] = c[i - 1] + a[i - 1];

            // update the safest index
            // what we want to do, is to be able to tell where the biggest
            // sum is so far on the fly
            // m[i] stores the index where the biggest sum is
            if (c[i - 1] < c[m[i - 1]])
            {
                m[i] = i - 1;
            }
            else
            {
                m[i] = m[i - 1];
            }

            // a(x, y) is the longest segment we found so far
            // what we want to test now, is whether there is a longer segment
            // than x, y, and that has to include i, since this is the only thing
            // that has changed
            int k = i - y + x - 1;
            while (k > 0)
            {
                if (c[i] - c[m[k]] >= sumL)
                {
                    k = m[k];
                }
                else
                {
                    break;
                }
                x = k + 1;
                y = i;
            }
        }

        if (x == 0)
        {
            return null;
        }

        return Range.between(x - 1, y - 1);
    }
}
