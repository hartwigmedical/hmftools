package com.hartwig.hmftools.redux.bqr;

import static com.hartwig.hmftools.redux.common.Constants.BQR_DUAL_AD;
import static com.hartwig.hmftools.redux.common.Constants.BQR_DUAL_AF_HIGH;
import static com.hartwig.hmftools.redux.common.Constants.BQR_DUAL_AF_LOW;
import static com.hartwig.hmftools.redux.common.Constants.BQR_NON_DUAL_AD;
import static com.hartwig.hmftools.redux.common.Constants.BQR_NON_DUAL_AF_HIGH;
import static com.hartwig.hmftools.redux.common.Constants.BQR_NON_DUAL_AF_LOW;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrReadType;

public class BaseQualityData
{
    public final byte[] TrinucleotideContext;

    private final Map<BqrReadType,List<AltQualityCount>> mAltQualityCountsMap;
    private boolean mHasIndel;

    public BaseQualityData(final byte[] trinucleotideContext)
    {
        TrinucleotideContext = trinucleotideContext;

        mHasIndel = false;
        mAltQualityCountsMap = Maps.newHashMap();
    }

    public byte ref() { return TrinucleotideContext[1]; }

    public void processReadBase(final BqrReadType readType, byte alt, byte quality, boolean posStrand)
    {
        List<AltQualityCount> altQualityCounts = mAltQualityCountsMap.get(readType);

        if(altQualityCounts == null)
        {
            altQualityCounts = Lists.newArrayList();
            mAltQualityCountsMap.put(readType, altQualityCounts);
        }

        for(AltQualityCount altQualityCount : altQualityCounts)
        {
            if(altQualityCount.Alt == alt && altQualityCount.Quality == quality)
            {
                altQualityCount.increment(posStrand);
                return;
            }
        }

        altQualityCounts.add(new AltQualityCount(alt, quality, posStrand));
    }

    public void setHasIndel() { mHasIndel = true; }
    public boolean hasIndel() { return mHasIndel; }

    public Map<BqrKey,Integer> formKeyCounts()
    {
        Map<BqrKey,Integer> keyCounts = Maps.newHashMap();

        // exclude any alt with too much support (regardless of quality)
        Map<Byte,Integer> standardAltCounts = Maps.newHashMap();
        Map<Byte,Integer> dualAltCounts = Maps.newHashMap();

        int standardTotalCount = 0;
        int dualTotalCount = 0;
        Set<Byte> allAlts = Sets.newHashSet();

        for(Map.Entry<BqrReadType,List<AltQualityCount>> entry : mAltQualityCountsMap.entrySet())
        {
            boolean isDual = entry.getKey().isHighQuality();

            for(AltQualityCount aqCount : entry.getValue())
            {
                int altTotalCount = aqCount.totalCount();

                if(isDual)
                    dualTotalCount += altTotalCount;
                else
                    standardTotalCount += altTotalCount;

                if(aqCount.Alt == ref())
                    continue;

                allAlts.add(aqCount.Alt);

                Map<Byte,Integer> altCountMap = isDual ? dualAltCounts : standardAltCounts;
                Integer altCount = altCountMap.get(aqCount.Alt);
                altCountMap.put(aqCount.Alt, altCount != null ? altCount + altTotalCount : altTotalCount);
            }
        }

        Set<Byte> failedAlts = Sets.newHashSet();

        for(Byte alt : allAlts)
        {
            boolean includeAlt = true;

            Integer standardAltCount = standardAltCounts.get(alt);

            if(standardAltCount != null)
            {
                double standardAltVaf = standardAltCount / (double)standardTotalCount;

                includeAlt = passesAltTest(standardAltCount, standardAltVaf, BQR_NON_DUAL_AF_LOW, BQR_NON_DUAL_AF_HIGH, BQR_NON_DUAL_AD);
            }

            if(includeAlt && dualAltCounts.containsKey(alt))
            {
                int dualAltCount = dualAltCounts.get(alt);
                double dualAltVaf = dualAltCount / (double)dualTotalCount;

                // for the dual condition it means: use a site if (AF<1% | AD<3) & AF <7.5%, or equivalently, AF<1% | (AD<3 & AF<7.5%)
                includeAlt = passesAltTest(dualAltCount, dualAltVaf, BQR_DUAL_AF_LOW, BQR_DUAL_AF_HIGH, BQR_DUAL_AD);
            }

            if(!includeAlt)
                failedAlts.add(alt);
        }

        for(Map.Entry<BqrReadType,List<AltQualityCount>> entry : mAltQualityCountsMap.entrySet())
        {
            BqrReadType readType = entry.getKey();

            for(AltQualityCount aqCount : entry.getValue())
            {
                if(failedAlts.contains(aqCount.Alt))
                    continue;

                if(aqCount.PosStrandCount > 0)
                {
                    keyCounts.put(new BqrKey(ref(), aqCount.Alt, TrinucleotideContext, aqCount.Quality, readType), aqCount.PosStrandCount);
                }

                if(aqCount.NegStrandCount > 0)
                {
                    byte[] tncReversed = Nucleotides.reverseComplementBases(TrinucleotideContext);
                    byte altReversed = Nucleotides.swapDnaBase(aqCount.Alt);
                    keyCounts.put(new BqrKey(tncReversed[1], altReversed, tncReversed, aqCount.Quality, readType), aqCount.NegStrandCount);
                }
            }
        }

        return keyCounts;
    }

    private static boolean passesAltTest(int altCount, double altVaf, double lowAfLimit, double highAfLimit, int adLimit)
    {
        return altVaf < lowAfLimit || (altVaf < highAfLimit && altCount <= adLimit);
    }

    public String toString()
    {
        return String.format("ref(%s) context(%s) readTypes(%d) alts(%d)",
                (char)ref(), new String(TrinucleotideContext), mAltQualityCountsMap.size(),
                mAltQualityCountsMap.values().stream().mapToInt(x -> x.size()).sum());
    }
}
