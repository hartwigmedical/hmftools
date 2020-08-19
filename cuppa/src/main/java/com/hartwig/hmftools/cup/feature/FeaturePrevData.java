package com.hartwig.hmftools.cup.feature;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

public class FeaturePrevData
{
    public final String CancerType;
    public final String Gene;
    public final FeatureType Type;
    public double Prevalence; // can be adjusted

    public FeaturePrevData(final String cancerType, final String gene, final FeatureType type, final double prevalence)
    {
        CancerType = cancerType;
        Gene = gene;
        Type = type;
        Prevalence = prevalence;
    }

    public static FeaturePrevData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 4)
            return null;

        return new FeaturePrevData(items[0], items[1], FeatureType.valueOf(items[2]), Double.parseDouble(items[3]));
    }

    public String toString() { return String.format("ct(%s) gene(%s) type(%s) prev(%.4f)", CancerType, Gene, Type, Prevalence); }

}
