package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

public class SampleDriverData
{
    public final String SampleId;
    public final String Gene;
    public final DriverType Type;

    public SampleDriverData(final String sampleId, final String gene, final DriverType type)
    {
        SampleId = sampleId;
        Gene = gene;
        Type = type;
    }

    public static SampleDriverData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 3)
            return null;

        return new SampleDriverData(items[0], items[1], DriverType.valueOf(items[2]));
    }

    public String toString() { return String.format("sample(%s) gene(%s) type(%s)", SampleId, Gene, Type); }
}
