package com.hartwig.hmftools.purple.drivers;

import com.hartwig.hmftools.common.region.BaseRegion;

public class AmpDelRegionFrequency
{
    public enum EventType
    {
        AMP,
        DEL
    }

    public static final String FLD_FREQUENCY = "Frequency";
    public static final String FLD_TYPE = "Type";

    public final BaseRegion Region;
    public int Frequency;
    public final EventType Type;

    public AmpDelRegionFrequency(final BaseRegion region, final EventType type, final int frequency)
    {
        Region = region;
        Type = type;
        Frequency = frequency;
    }
}
