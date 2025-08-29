package com.hartwig.hmftools.common.segmentation;

/**
 * Calculates a segment penalty based on an array of values.
 */
public class Gamma
{
    private final double segmentPenalty;

    public Gamma(double[] y, double gamma, boolean normalise)
    {
        if(y.length == 0)
        {
            throw new IllegalArgumentException("Input array must not be empty");
        }

        if(normalise)
        {
            double[] runningMedians = new Runmed().runmed(y, filterWidth(y.length));

            double[] diffs = new double[y.length];
            for(int i = 0; i < y.length; i++)
            {
                diffs[i] = y[i] - runningMedians[i];
            }
            double SD = Stats.medianAbsoluteDeviation(diffs);
            segmentPenalty = SD * SD * gamma;
        }
        else
        {
            segmentPenalty = gamma;
        }
    }

    public Gamma(double[] y, double gamma)
    {
        this(y, gamma, false);
    }

    public Gamma(double[] y)
    {
        this(y, 50.0, false);
    }

    public double getSegmentPenalty()
    {
        return segmentPenalty;
    }

    public static int filterWidth(int n)
    {
        return n >= 51 ? 51 : (n % 2 == 0 ? n - 1 : n);
    }
}