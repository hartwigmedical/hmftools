package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class RawCounters
{
    public int RawDepth;
    public int RawAltSupport;
    public int RawRefSupport;
    public int RawAltBaseQualityTotal;
    public int RawRefBaseQualityTotal;
    public int RawContextAltBaseQualityTotal;

    // TEMP: diffs between the read context match result and raw context
    public int RawVsMatchRefDiffs;
    public int RawVsMatchAltDiffs;

    public RawCounters()
    {
        RawDepth = 0;
        RawAltSupport = 0;
        RawRefSupport = 0;
        RawAltBaseQualityTotal = 0;
        RawRefBaseQualityTotal = 0;
        RawContextAltBaseQualityTotal = 0;
        RawVsMatchRefDiffs = 0;
        RawVsMatchAltDiffs = 0;
    }

    public void registerRawSupport(final RawContext rawContext)
    {
        if(rawContext.AltSupport)
        {
            ++RawAltSupport;
            RawAltBaseQualityTotal += rawContext.BaseQuality;
        }
        else if(rawContext.RefSupport)
        {
            ++RawRefSupport;
            RawRefBaseQualityTotal += rawContext.BaseQuality;
        }
    }

    public String toString()
    {
        return format("depth(%d) support(ref=%d alt=%d) qual(ref=%d alt=%d raw=%d)",
                RawDepth, RawRefSupport, RawAltSupport, RawRefBaseQualityTotal, RawAltBaseQualityTotal, RawRefBaseQualityTotal);
    }
}
