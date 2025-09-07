package com.hartwig.hmftools.common.segmentation;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

class Stats
{
    static DescriptiveStatistics ds(double[] array)
    {
        DescriptiveStatistics ds = new DescriptiveStatistics();
        for(double value : array)
        {
            ds.addValue(value);
        }
        return ds;
    }

    static double mean(double[] array)
    {
        return ds(array).getMean();
    }

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
        DescriptiveStatistics ds = ds(array);
        double median = ds.getPercentile(50.0);
        DescriptiveStatistics diffsDS = new DescriptiveStatistics();
        for(double value : array)
        {
            diffsDS.addValue(Math.abs(value - median));
        }
        return diffsDS.getPercentile(50.0) * 1.4826;
    }
}