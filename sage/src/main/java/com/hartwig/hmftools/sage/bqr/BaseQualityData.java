package com.hartwig.hmftools.sage.bqr;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class BaseQualityData
{
    public final int Position;
    public final byte Ref;
    public final byte[] TrinucleotideContext;

    private final List<AltQualityCount> mAltQualityCounts;
    private boolean mHasIndel;

    public BaseQualityData(final int position, final byte ref, final byte[] trinucleotideContext)
    {
        Position = position;
        Ref = ref;
        TrinucleotideContext = trinucleotideContext;

        mHasIndel = false;
        mAltQualityCounts = Lists.newArrayList();
    }

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

    public Map<BqrKey,Integer> formKeyCounts(int maxAltCount, double maxAltPerc)
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

        for(AltQualityCount aqCount : mAltQualityCounts)
        {
            if(altCounts.containsKey(aqCount.Alt))
            {
                int altCount = altCounts.get(aqCount.Alt);
                double altVaf = altCount / (double)totalCount;
                if(altVaf > maxAltPerc && altCount > maxAltCount)
                    continue;
            }

            keyCounts.put(new BqrKey(Ref, aqCount.Alt, TrinucleotideContext, aqCount.Quality), aqCount.Count);
        }

        return keyCounts;
    }

    public String toString()
    {
        return String.format("%d: ref(%s) context(%s) alts(%d)",
                Position, (char)Ref, new String(TrinucleotideContext), mAltQualityCounts.size());
    }

    private class AltQualityCount
    {
        public final byte Alt;
        public final byte Quality;
        public int Count;

        public AltQualityCount(final byte alt, final byte quality)
        {
            Alt = alt;
            Quality = quality;
            Count = 1;
        }

        public String toString() { return String.format("alt(%s) qual(%d) count(%d)", (char)Alt, (int)Quality, Count); }
    }
}
