package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.min;

public class Percentiles
{
    public static final double INVALID_VALUE = -1;
    public static final int PERCENTILE_COUNT = 101;

    public static double getPercentile(final double[] percentileValues, double value)
    {
        // find the position of a value within a set of percentile values
        if(percentileValues == null || percentileValues.length != PERCENTILE_COUNT)
            return INVALID_VALUE;

        if(value < percentileValues[0])
            return 0;
        else if(value > percentileValues[percentileValues.length - 1])
            return (percentileValues.length - 1) * 0.01;

        for(int i = 0; i < percentileValues.length - 1; ++i)
        {
            if(value >= percentileValues[i] && value <= percentileValues[i + 1])
            {
                if(percentileValues[i + 1] == percentileValues[i])
                    return i * 0.01;

                double upperFactor = (value - percentileValues[i]) / (percentileValues[i + 1] - percentileValues[i]);
                return (upperFactor * (i + 1) + (1 - upperFactor) * i) * 0.01;
            }
        }

        return INVALID_VALUE;
    }

    public static double[] buildPercentiles(final double[] values)
    {
        double[] percentiles = new double[PERCENTILE_COUNT];
        calcPercentileValues(values, percentiles);
        return percentiles;
    }

    public static void calcPercentileValues(final double[] values, final double[] percentileValues)
    {
        // assumes that values are sorted in ascending order
        int sampleCount = values.length;

        // populate the upper and lower bounds
        double percSlots = percentileValues.length;

        double samplesPerPercentile = sampleCount/percSlots;

        for(int i = 0; i < percentileValues.length; ++i)
        {
            double lowerIndex = i * samplesPerPercentile;
            double upperIndex = lowerIndex + samplesPerPercentile * 0.9999;

            int lowerBound = (int)floor(lowerIndex);
            int upperBound = (int)ceil(upperIndex) - 1;
            upperBound = min(upperBound, sampleCount);

            if(lowerBound == upperBound)
            {
                percentileValues[i] = values[lowerBound];
                continue;
            }

            double tpmTotal = 0;
            double sampleTotal = 0;

            for(int s = lowerBound; s <= upperBound; ++s)
            {
                double tpm = values[s];

                double fractionOfTpm;

                if(s == lowerBound)
                {
                    fractionOfTpm = 1 - (lowerIndex - lowerBound);
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else if(s == upperBound)
                {
                    fractionOfTpm = upperIndex - upperBound;
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else
                {
                    ++sampleTotal;
                    tpmTotal += tpm;
                }
            }

            percentileValues[i] = tpmTotal / sampleTotal;
        }
    }

}
