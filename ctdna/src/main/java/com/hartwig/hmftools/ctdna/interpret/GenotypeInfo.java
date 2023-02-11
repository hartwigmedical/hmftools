package com.hartwig.hmftools.ctdna.interpret;

import static java.lang.String.format;

public class GenotypeInfo
{
    public final int Index;
    public final String Sample;

    public GenotypeInfo(final int index, final String sample)
    {
        Index = index;
        Sample = sample;
    }

    public String toString() { return format("%d:%s", Index, Sample); }
}
