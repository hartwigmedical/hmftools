package com.hartwig.hmftools.common.segmentation.copynumber;

/**
 * Java implementation of R's runmed function.
 * The R function has options for dealing with the endpoints but in
 * copynumber the only option used is "median".
 * <p>
 * <a href="https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/runmed.R">...</a>
 */
public class Runmed
{
    public double[] runmed(double[] data, int k, boolean smooth)
    {
        if(k < 1)
        {
            throw new IllegalArgumentException("Window size k must be >= 1");
        }
        if(k % 2 != 1)
        {
            throw new IllegalArgumentException("Window size k must be odd");
        }
        if(data.length < 2)
        {
            return data.clone();
        }

        double[] result = new WindowedMedian(data, k).getMedians();

        return smooth ? smoothEnds(result, k) : result;
    }

    public double[] runmed(double[] data, int k)
    {
        return runmed(data, k, true);
    }

    public static double[] smoothEnds(double[] y, int k)
    {
        if(k < 0 || k % 2 == 0)
        {
            throw new IllegalArgumentException("bandwidth 'k' must be >= 1 and odd!");
        }

        int halfK = k / 2;
        if(halfK < 1)
        {
            return y.clone(); // Nothing to do if k < 3
        }

        // Clone the input array to avoid modifying it
        int n = y.length;
        double[] sm = y.clone();

        // Apply smoothing for k >= 3
        if(halfK >= 2)
        {
            // Apply med3 for the second and second-to-last points
            sm[1] = med3(y[0], y[1], y[2]);
            sm[n - 2] = med3(y[n - 1], y[n - 2], y[n - 3]);

            // For larger k values, use medianOdd for the border points
            if(halfK >= 3)
            {
                for(int i = 3; i <= halfK; i++)
                {
                    int j = 2 * i - 1;
                    sm[i - 1] = medianOdd(sliceArray(y, 0, j), j); // Left border
                    sm[n - i] = medianOdd(sliceArray(y, n - j, n), j); // Right border
                }
            }
        }

        // Apply Tukey's end-point rule for the very first and last points
        // For the first point: median of (original value, second point value, extrapolated value)
        // The extrapolated value is calculated as 3*sm[1] - 2*sm[2]
        sm[0] = med3(y[0], sm[1], 3 * sm[1] - 2 * sm[2]);
        sm[n - 1] = med3(y[n - 1], sm[n - 2], 3 * sm[n - 2] - 2 * sm[n - 3]);

        return sm;
    }

    public static double med3(double a, double b, double c)
    {
        if(a < b)
        { //...a...b...
            if(c < b)
            {
                return Math.max(a, c);
            }
            else
            {
                return b;
            }
        }
        else
        { //...b...a...
            if(c < a)
            {
                return Math.max(c, b);
            }
            else
            {
                return a;
            }
        }
    }

    /**
     * Computes the median of an array of doubles when n is odd.
     * This is slightly more efficient than a general median calculation.
     */
    public static double medianOdd(double[] x, int n)
    {
        if(n % 2 != 1)
        {
            throw new IllegalArgumentException("n must be odd for medianOdd");
        }
        double[] sorted = x.clone();
        java.util.Arrays.sort(sorted);
        int half = (n + 1) / 2;
        return sorted[half - 1];
    }

    private static double[] sliceArray(double[] array, int start, int end)
    {
        double[] result = new double[end - start];
        System.arraycopy(array, start, result, 0, end - start);
        return result;
    }
}