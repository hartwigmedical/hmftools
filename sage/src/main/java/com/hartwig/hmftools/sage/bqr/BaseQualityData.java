package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageConstants.BQR_DUAL_AD;
import static com.hartwig.hmftools.sage.SageConstants.BQR_DUAL_AF_HIGH;
import static com.hartwig.hmftools.sage.SageConstants.BQR_DUAL_AF_LOW;
import static com.hartwig.hmftools.sage.SageConstants.BQR_NON_DUAL_AD;
import static com.hartwig.hmftools.sage.SageConstants.BQR_NON_DUAL_AF_HIGH;
import static com.hartwig.hmftools.sage.SageConstants.BQR_NON_DUAL_AF_LOW;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;

public class BaseQualityData
{
    public final byte Ref;
    public final byte[] TrinucleotideContext;
    public final BqrReadType ReadType;

    private final List<AltQualityCount> mAltQualityCounts;
    private boolean mHasIndel;

    public BaseQualityData(final byte ref, final byte[] trinucleotideContext, final BqrReadType readType)
    {
        Ref = ref;
        TrinucleotideContext = trinucleotideContext;
        ReadType = readType;

        mHasIndel = false;
        mAltQualityCounts = Lists.newArrayList();
    }

    public List<AltQualityCount> altQualityCounts() { return mAltQualityCounts; }

    public void processReadBase(byte alt, byte quality)
    {
        for(AltQualityCount altQualityCount : mAltQualityCounts)
        {
            if(altQualityCount.Alt == alt && altQualityCount.Quality == quality)
            {
                ++altQualityCount.Count;
                return;
            }
        }

        mAltQualityCounts.add(new AltQualityCount(alt, quality));
    }

    public void setHasIndel() { mHasIndel = true; }
    public boolean hasIndel() { return mHasIndel; }

    public Map<BqrKey,Integer> formKeyCounts()
    {
        Map<BqrKey,Integer> keyCounts = Maps.newHashMap();

        // exclude any alt with too much support (regardless of quality)
        Map<Byte,Integer> altCounts = Maps.newHashMap();

        int totalCount = 0;
        for(AltQualityCount aqCount : mAltQualityCounts)
        {
            totalCount += aqCount.Count;

            if(aqCount.Alt == Ref)
                continue;

            Integer altCount = altCounts.get(aqCount.Alt);
            altCounts.put(aqCount.Alt, altCount != null ? altCount + aqCount.Count : aqCount.Count);
        }

        double lowAfLimit = ReadType.isHighQuality() ? BQR_DUAL_AF_LOW : BQR_NON_DUAL_AF_LOW;
        double highAfLimit = ReadType.isHighQuality() ? BQR_DUAL_AF_HIGH : BQR_NON_DUAL_AF_HIGH;
        int adLimit = ReadType.isHighQuality() ? BQR_DUAL_AD : BQR_NON_DUAL_AD;

        for(AltQualityCount aqCount : mAltQualityCounts)
        {
            if(altCounts.containsKey(aqCount.Alt))
            {
                int altCount = altCounts.get(aqCount.Alt);
                double altVaf = altCount / (double)totalCount;

                // for the dual condition it means: use a site if (AF<1% | AD<3) & AF <7.5%, or equivalently, AF<1% | (AD<3 & AF<7.5%)
                boolean includeAlt = altVaf < lowAfLimit || (altVaf < highAfLimit && altCount <= adLimit);

                if(!includeAlt)
                    continue;
            }

            keyCounts.put(new BqrKey(Ref, aqCount.Alt, TrinucleotideContext, aqCount.Quality, ReadType), aqCount.Count);
        }

        return keyCounts;
    }

    public String toString()
    {
        return String.format("ref(%s) context(%s) readType(%s) alts(%d)",
                (char)Ref, new String(TrinucleotideContext), ReadType, mAltQualityCounts.size());
    }

}
