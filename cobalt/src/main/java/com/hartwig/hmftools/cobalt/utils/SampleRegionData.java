package com.hartwig.hmftools.cobalt.utils;

import static java.lang.String.format;

public class SampleRegionData
{
    public final int ReadCount;
    public final double GcRatioPanel;
    public final double GcRatioWgs;

    public SampleRegionData(final int readCount, final double gcRatioPanel, final double gcRatioWgs)
    {
        ReadCount = readCount;
        GcRatioPanel = gcRatioPanel;
        GcRatioWgs = gcRatioWgs;
    }

    public String toString() { return format("reads(%d) ratio(panel=%.3f wgs=%f)", ReadCount, GcRatioPanel, GcRatioWgs); }
}
