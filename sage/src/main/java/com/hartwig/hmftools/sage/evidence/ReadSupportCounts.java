package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.variant.VariantReadSupport;

public class ReadSupportCounts
{
    public int Full;
    public int PartialCore;
    public int Core;
    public int Realigned;
    public int Ref;
    public int Total;

    public ReadSupportCounts()
    {
        Full = 0;
        PartialCore = 0;
        Core = 0;
        Realigned = 0;
        Ref = 0;
        Total = 0;
    }

    public void addSupport(@Nullable final VariantReadSupport support, final int count)
    {
        Total += count;

        if(support != null)
        {
            switch(support)
            {
                case FULL:
                    Full += count;
                    break;

                case PARTIAL_CORE:
                    PartialCore += count;
                    break;

                case CORE:
                    Core += count;
                    break;

                case REALIGNED:
                    Realigned += count;
                    break;

                case REF:
                    Ref += count;
                    break;
            }
        }
    }

    public int altSupport() { return Full + PartialCore + Core + Realigned; }
    public int strongSupport() { return Full + PartialCore + Realigned; }

    public void applyRatio(double ratio)
    {
        if(ratio == 1)
            return;

        Full = (int)round(Full * ratio);
        PartialCore = (int)round(PartialCore * ratio);
        Core = (int)round(Core * ratio);
        Realigned = (int)round(Realigned * ratio);
        Ref = (int)round(Ref * ratio);
        Total = (int)round(Total * ratio);
    }

    public int[] toArray()
    {
        return new int[] { Full, PartialCore, Core, Realigned, Ref, Total };
    }
}
