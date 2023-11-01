package com.hartwig.hmftools.sage.evidence;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.variant.VariantReadSupport;

public class ReadSupportCounts
{
    public int Full;
    public int Partial;
    public int Core;
    public int Realigned;
    public int OtherAlt;
    public int Ref;
    public int Total;

    public ReadSupportCounts()
    {
        Full = 0;
        Partial = 0;
        Core = 0;
        Realigned = 0;
        OtherAlt = 0;
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

                case PARTIAL:
                    Partial += count;
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

                case SIMPLE_ALT:
                    OtherAlt += count;
                    break;
            }
        }
    }

    public int altSupport() { return Full + Partial + Core + OtherAlt + Realigned; }

    public int[] toArray()
    {
        return new int[] { Full, Partial, Core, Realigned, OtherAlt, Ref, Total };
    }

}
