package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

public class DriverPrevData
{
    public final String CancerType;
    public final String Gene;
    public final DriverType Type;
    public final double Prevalence;

    public DriverPrevData(final String cancerType, final String gene, final DriverType type, final double prevalence)
    {
        CancerType = cancerType;
        Gene = gene;
        Type = type;
        Prevalence = prevalence;
    }

    public static DriverPrevData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 4)
            return null;

        return new DriverPrevData(items[0], items[1], DriverType.valueOf(items[2]), Double.parseDouble(items[3]));
    }

    public String toString() { return String.format("ct(%s) gene(%s) type(%s) prev(%.4f)", CancerType, Gene, Type, Prevalence); }

}
