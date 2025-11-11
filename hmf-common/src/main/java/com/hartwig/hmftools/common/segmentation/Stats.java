package com.hartwig.hmftools.common.segmentation;

import com.hartwig.hmftools.common.utils.Doubles;

class Stats
{
    /**
     * Computes the median absolute deviation of an array of doubles.
     * Copied from the R code: <a href="https://search.r-project.org/R/refmans/stats/html/mad.html">...</a>
     * <p>
     * Arguments
     * x a numeric vector.
     * center Optionally, the centre: defaults to the median.
     * constant scale factor.
     * na.rm if TRUE then NA values are stripped from x before computation takes place.
     * low if TRUE, compute the 'lo-median', i.e., for even sample size, do not average the two middle values, but take the smaller one.
     * high if TRUE, compute the 'hi-median', i.e., take the larger of the two middle values for even sample size.
     * <p>
     * Details
     * The actual value calculated is constant * cMedian(abs(x - center)) with the default value of center being median(x), and cMedian being the usual, the 'low' or 'high' median.
     * <p>
     * The default constant = 1.4826 (approximately...)
     */
    static double medianAbsoluteDeviation(double[] array)
    {
        double median = Doubles.median(array);
        double[] absoluteDeviations = new double[array.length];
        for(int i = 0; i < array.length; i++)
        {
            absoluteDeviations[i] = Math.abs(array[i] - median);
        }
        return Doubles.median(absoluteDeviations) * 1.4826;
    }
}
