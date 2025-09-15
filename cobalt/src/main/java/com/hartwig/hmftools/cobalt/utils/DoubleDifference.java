package com.hartwig.hmftools.cobalt.utils;

public class DoubleDifference
{
    public final Double difference;
    public final boolean hasDifference;

    public DoubleDifference(Double original, Double comparison, Double epsilon)
    {
        double diff = original - comparison;
        if(Math.abs(diff) < epsilon)
        {
            difference = 0.0;
            hasDifference = false;
        }
        else
        {
            difference = diff;
            hasDifference = true;
        }
    }
}
