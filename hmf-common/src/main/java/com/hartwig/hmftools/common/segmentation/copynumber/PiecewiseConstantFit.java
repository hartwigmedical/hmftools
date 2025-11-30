package com.hartwig.hmftools.common.segmentation.copynumber;

import java.util.Arrays;

public record PiecewiseConstantFit(int[] lengths, int[] startPositions, double[] means)
{
    @Override
    public boolean equals(Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        PiecewiseConstantFit pcf = (PiecewiseConstantFit) o;

        if(!Arrays.equals(lengths, pcf.lengths))
        {
            return false;
        }
        if(!Arrays.equals(startPositions, pcf.startPositions))
        {
            return false;
        }
        return Arrays.equals(means, pcf.means);
    }

    @Override
    public int hashCode()
    {
        int result = Arrays.hashCode(lengths);
        result = 31 * result + Arrays.hashCode(startPositions);
        result = 31 * result + Arrays.hashCode(means);
        return result;
    }

    @Override
    public String toString()
    {
        return "PiecewiseConstantFit{" +
                "lengths=" + Arrays.toString(lengths) +
                ", startPositions=" + Arrays.toString(startPositions) +
                ", means=" + Arrays.toString(means) +
                '}';
    }
}