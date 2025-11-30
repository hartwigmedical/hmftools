package com.hartwig.hmftools.common.segmentation;

import com.google.common.base.Preconditions;

public record DataForSegmentation(double[] valuesForSegmentation, double[] rawValues)
{
    public DataForSegmentation
    {
        Preconditions.checkArgument(valuesForSegmentation.length == rawValues.length);
    }

    public int count()
    {
        return valuesForSegmentation.length;
    }

    public boolean isEmpty()
    {
        return count() == 0;
    }
}
