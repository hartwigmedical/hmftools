package com.hartwig.hmftools.common.sigs;

public class Percentiles
{
    public static final double INVALID_VALUE = -1;
    private static final int PERCENTILE_COUNT = 101;

    public static double getPercentile(final double[] percentileValues, double value)
    {
        if(percentileValues == null || percentileValues.length != PERCENTILE_COUNT)
            return INVALID_VALUE;

        if(value < percentileValues[0])
            return 0;
        else if(value > percentileValues[percentileValues.length - 1])
            return percentileValues.length - 1;

        for(int i = 0; i < percentileValues.length - 1; ++i)
        {
            if(value >= percentileValues[i] && value <= percentileValues[i + 1])
            {
                if(percentileValues[i + 1] == percentileValues[i])
                    return i;

                double upperFactor = (value - percentileValues[i]) / (percentileValues[i + 1] - percentileValues[i]);
                return upperFactor * (i + 1) + (1 - upperFactor) * i;
            }
        }

        return INVALID_VALUE;
    }

    public double[] buildPercentiles(final double[] values)
    {
        double[] percentiles = new double[PERCENTILE_COUNT];
        return percentiles;
    }

}
