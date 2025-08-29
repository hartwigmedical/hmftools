package com.hartwig.hmftools.common.segmentation;

import java.util.Arrays;

public class PCF
{
    private final int[] lengths;
    private final int[] startPositions;
    private final double[] means;

    public PCF(int[] lengths, int[] startPositions, double[] means)
    {
        this.lengths = lengths;
        this.startPositions = startPositions;
        this.means = means;
    }

    public int[] getLengths()
    {
        return lengths;
    }

    public int[] getStartPositions()
    {
        return startPositions;
    }

    public double[] getMeans()
    {
        return means;
    }

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

        PCF pcf = (PCF) o;

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
}