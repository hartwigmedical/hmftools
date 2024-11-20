package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.bqr.BqrRegionReader.extractReadType;
import static java.lang.Math.round;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import htsjdk.samtools.SAMRecord;

public class ReadSupportCounts
{
    public int Full;
    public int PartialCore;
    public int Core;
    public int Realigned;
    public int Ref;
    public int Total;
    public int StrongSimplexSupport;
    public int FullDuplex;

    public ReadSupportCounts()
    {
        Full = 0;
        PartialCore = 0;
        Core = 0;
        Realigned = 0;
        Ref = 0;
        Total = 0;
        StrongSimplexSupport = 0;
        FullDuplex = 0;
    }

    public void addSupport(SAMRecord record, int readIndex, @Nullable final VariantReadSupport support, final int count)
    {
        Total += count;
        if(support == VariantReadSupport.FULL || support == VariantReadSupport.PARTIAL_CORE || support == VariantReadSupport.REALIGNED)
        {
            BqrReadType readType = extractReadType(record, SequencingType.SBX, record.getBaseQualities()[readIndex]);
            if(support == VariantReadSupport.FULL && readType == BqrReadType.DUAL)
                FullDuplex += count;
            if(readType != BqrReadType.DUAL)
                StrongSimplexSupport += count;
        }

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
    public int strongSimplexSupport() { return StrongSimplexSupport; }

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
