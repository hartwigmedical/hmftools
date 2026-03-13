package com.hartwig.hmftools.qsee.common;

import java.util.Arrays;

public class LinearInterpolation
{
    public static double interpolate1D(double lowerValue, double upperValue, double fraction)
    {
        return lowerValue + fraction * (upperValue - lowerValue);
    }

    public static double interpolate2D(double inputXValue, double[] XValues, double[] YValues)
    {
        if(inputXValue < XValues[0])
            return Double.NEGATIVE_INFINITY;

        if(inputXValue > XValues[XValues.length - 1])
            return Double.POSITIVE_INFINITY;

        int matchIndex = Arrays.binarySearch(XValues, inputXValue);
        if(matchIndex >= 0)
        {
            return YValues[matchIndex];
        }

        int upperIndex = -matchIndex - 1; // When no exact match is found with binary search, a negative match index is returned encoding the insertion point
        int lowerIndex =  upperIndex - 1;

        double lowerXValue = XValues[lowerIndex];
        double upperXValue = XValues[upperIndex];
        double lowerYValue = YValues[lowerIndex];
        double upperYValue = YValues[upperIndex];

        if(lowerXValue == inputXValue)
        {
            return YValues[lowerIndex];
        }

        double distanceToAdd = inputXValue - lowerXValue;
        double distanceTotal = upperXValue - lowerXValue;
        double fraction = distanceToAdd / distanceTotal;
        return interpolate1D(lowerYValue, upperYValue, fraction);
    }
}
