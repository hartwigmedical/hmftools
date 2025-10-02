package com.hartwig.hmftools.cobalt.utils;

public class DoubleDifference
{
    public final Double difference;
    public final boolean hasDifference;

    public DoubleDifference(Double original, Double comparison, Double epsilon)
    {
        if (original == 0 && comparison == 0)
        {
            difference = 0.0;
            hasDifference = false;
        }
        else
        {
            double diff = ((original - comparison) * 2) / (Math.abs(original) + Math.abs(comparison));
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
}
