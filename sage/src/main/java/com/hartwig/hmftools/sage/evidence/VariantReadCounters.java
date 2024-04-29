package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class VariantReadCounters
{
    public int Depth;
    public int RecalibratedBaseQualityTotal;
    public double SupportAltBaseQualityTotal;
    public long MapQualityTotal;
    public long AltMapQualityTotal;

    // qual totals for raw (ie non-model) base quality
    public int AltRawBaseQualityTotal;
    public int RefRawBaseQualityTotal;

    public VariantReadCounters()
    {
        RecalibratedBaseQualityTotal = 0;
        SupportAltBaseQualityTotal = 0;
        MapQualityTotal = 0;
        AltMapQualityTotal = 0;

        Depth = 0;
        AltRawBaseQualityTotal = 0;
        RefRawBaseQualityTotal = 0;
    }

    /*
    public void registerRawSupport(final double rawBaseQuality)
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
    */

    public String toString()
    {
        return format("depth(%d) qual(ref=%d alt=%d)",
                Depth, RefRawBaseQualityTotal, AltRawBaseQualityTotal);
    }
}
