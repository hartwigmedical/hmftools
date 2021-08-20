package com.hartwig.hmftools.common.stats;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.capValue;

public final class Percentiles
{
    public static final double INVALID_VALUE = -1;
    public static final int PERCENTILE_COUNT = 101;

    public static double getPercentile(final double[] percentileValues, double value)
    {
        return getPercentile(percentileValues, value, false, 0);
    }

    public static double getPercentile(final double[] percentileValues, double value, boolean useMaxMultiple, double undefinedMaxMultiple)
    {
        // find the position of a value within a set of percentile values
        if(percentileValues == null || percentileValues.length != PERCENTILE_COUNT)
            return INVALID_VALUE;

        if(value < percentileValues[0])
        {
            // report values below the lowest value as a negative inverse fraction
            if(percentileValues[0] == 0 || !useMaxMultiple)
                return 0;
            else if(value == 0)
                return -100;
            else
                return -percentileValues[0] / value;
        }
        else if(value > percentileValues[percentileValues.length - 1])
        {
            // all values in the distribution are zero, so the multiple is undefined
            if(useMaxMultiple)
            {
                double maxValue = percentileValues[percentileValues.length - 1];
                return maxValue > 0 ? value / maxValue : undefinedMaxMultiple;
            }

            return (percentileValues.length - 1) * 0.01;
        }

        int sameValueStart = -1;
        int sameValueEnd = 0;
        for(int i = 0; i < percentileValues.length - 1; ++i)
        {
            // where successive percentile values are the same (eg zeros for the first X percentiles), report the median percentile
            if(value == percentileValues[i] && percentileValues[i] == percentileValues[i + 1])
            {
                if(sameValueStart == -1)
                    sameValueStart =  i;

                sameValueEnd = i + 1;

                if(i < percentileValues.length - 2)
                    continue;
            }

            if(value >= percentileValues[i] && value <= percentileValues[i + 1])
            {
                if(value > percentileValues[i] && value == percentileValues[i + 1] && i < percentileValues.length - 2)
                    continue;

                if(value == percentileValues[i] && sameValueStart >= 0 && sameValueStart < sameValueEnd)
                    return (0.5 * sameValueStart + 0.5 * sameValueEnd) * 0.01;

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
        // need to handle either #values being less than #percentiles (or slots if not using percentiles) and vice versa
        if(values == null || percentileValues == null)
            return;

        if(values.length == 1)
        {
            for(int i = 0; i < percentileValues.length; ++i)
            {
                percentileValues[i] = values[0];
            }

            return;
        }

        int valueCount = values.length;
        int slotCount = percentileValues.length;
        double valuesPerSlot = (valueCount - 1) / (double) (slotCount - 1);

        percentileValues[0] = values[0];

        for(int i = 1; i < percentileValues.length; ++i)
        {
            double valueIndex = i * valuesPerSlot;
            int valueLowerIndex = (int) floor(valueIndex);
            int valueUpperIndex = min((int) ceil(valueIndex), valueCount - 1);

            double valueLower = values[valueLowerIndex];
            double valueUpper = values[valueUpperIndex];
            double lowerFraction = capValue(1 - (valueIndex - valueLowerIndex), 0, 1);

            percentileValues[i] = lowerFraction * valueLower + (1 - lowerFraction) * valueUpper;
        }
    }

}
