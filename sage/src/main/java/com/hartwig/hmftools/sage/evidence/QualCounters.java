package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class QualCounters
{
    public int BaseQualityTotal;
    public int AltBaseQualityTotal;
    public long MapQualityTotal;
    public long AltMapQualityTotal;

    public QualCounters()
    {
        BaseQualityTotal = 0;
        AltBaseQualityTotal = 0;
        MapQualityTotal = 0;
        AltMapQualityTotal = 0;
    }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                BaseQualityTotal, AltBaseQualityTotal, MapQualityTotal, AltMapQualityTotal);
    }
}
