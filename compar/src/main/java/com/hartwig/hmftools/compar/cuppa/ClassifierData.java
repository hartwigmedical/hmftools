package com.hartwig.hmftools.compar.cuppa;

public class ClassifierData
{
    public final String DataType;
    public final String TopRefCancerType;
    public final double TopRefCancerValue;

    public ClassifierData(final String dataType, final String topRefCancerType, final double topRefCancerValue)
    {
        DataType = dataType;
        TopRefCancerType = topRefCancerType;
        TopRefCancerValue = topRefCancerValue;
    }
}
