package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_BOOSTED_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_T0;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_TP;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_0_BOOST;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractTpValues;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.coreHomopolymerLengths;
import static com.hartwig.hmftools.sage.vcf.ReadContextVcfInfo.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class UltimaVariantData
{
    private final VariantReadContext mReadContext;
    private final List<Integer> mHomopolymerLengths;
    private final List<Double> mHomopolymerAvgQuals;
    private final List<Double> mT0AvgQuals;
    private int mHomopolymerAvgQualsCount;

    public UltimaVariantData(final VariantReadContext readContext)
    {
        mReadContext = readContext;
        mHomopolymerLengths = coreHomopolymerLengths(readContext);

        mHomopolymerAvgQuals = Lists.newArrayList();
        mT0AvgQuals = Lists.newArrayList();

        while(mHomopolymerAvgQuals.size() < mHomopolymerLengths.size())
        {
            mHomopolymerAvgQuals.add(0.0);
            mT0AvgQuals.add(0.0);
        }
    }

    public List<Integer> homopolymerLengths() { return mHomopolymerLengths; }
    public List<Double> homopolymerAvgQuals() { return mHomopolymerAvgQuals; }
    public List<Double> t0AvgQuals() { return mT0AvgQuals; }

    public void addReadSupportInfo(final SAMRecord record, final int readVarIndex)
    {
        String recordCore = record.getReadString()
                .substring(readVarIndex - mReadContext.leftCoreLength(), readVarIndex + mReadContext.rightCoreLength() + 1);

        if(recordCore.equals(mReadContext.coreStr()))
        {
            registerHomopolymerQuals(record, readVarIndex);
            registerT0Quals(record, readVarIndex);
        }
    }

    private void registerHomopolymerQuals(final SAMRecord record, int readVarIndex)
    {
        byte[] baseQuals = record.getBaseQualities();
        byte[] tpValues = extractTpValues(record);
        int readIndex = readVarIndex - mReadContext.leftCoreLength();
        boolean firstHomopolymer = true;
        List<Integer> homopolyerQuals = Lists.newArrayList();

        for(int len : mHomopolymerLengths)
        {
            int lookupIndex = firstHomopolymer ? readIndex + len - 1 : readIndex;
            int tpValue = tpValues[lookupIndex];

            int homopolymerQual;
            if(tpValue == 0)
            {
                homopolymerQual = ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST;
            }
            else if(len == 1)
            {
                homopolymerQual = baseQuals[lookupIndex];
            }
            else
            {
                homopolymerQual = baseQuals[lookupIndex] - 3;
            }

            homopolyerQuals.add(homopolymerQual);

            firstHomopolymer = false;
            readIndex += len;
        }

        for(int i = 0; i < homopolyerQuals.size(); i++)
        {
            double avgQual = (mHomopolymerAvgQuals.get(i) * mHomopolymerAvgQualsCount + homopolyerQuals.get(i)) / (mHomopolymerAvgQualsCount + 1);
            mHomopolymerAvgQuals.set(i, avgQual);
        }
    }

    private void registerT0Quals(final SAMRecord record, int readVarIndex)
    {
        byte[] t0Values = extractT0Values(record);
        int readIndex = readVarIndex - mReadContext.leftCoreLength();
        boolean firstHomopolymer = true;
        List<Integer> t0Quals = Lists.newArrayList();
        for(int len : mHomopolymerLengths)
        {
            int lookupIndex = firstHomopolymer ? readIndex + len - 1 : readIndex;
            int t0Value = t0Values[lookupIndex];

            int t0Qual;
            if(t0Value == ULTIMA_BOOSTED_QUAL)
            {
                t0Qual = ULTIMA_MAX_QUAL_T0;
            }
            else
            {
                t0Qual = t0Value;
            }

            t0Quals.add(t0Qual);

            firstHomopolymer = false;
            readIndex += len;
        }

        for(int i = 0; i < t0Quals.size(); i++)
        {
            double avgQual = (mT0AvgQuals.get(i) * mHomopolymerAvgQualsCount + t0Quals.get(i)) / (mHomopolymerAvgQualsCount + 1);
            mT0AvgQuals.set(i, avgQual);
        }

        mHomopolymerAvgQualsCount++;
    }

    public String coreHomopolymerInfo()
    {
        if(mHomopolymerLengths == null)
            return null;

        String lengthsStr = mHomopolymerLengths.stream()
                .map(String::valueOf)
                .collect(Collectors.joining(ITEM_DELIM));

        String qualsStr = mHomopolymerAvgQuals.stream()
                .map(x -> format("%.3f", x))
                .collect(Collectors.joining(ITEM_DELIM));

        return lengthsStr + ITEM_DELIM + qualsStr;
    }

    public String coreT0Info()
    {
        if(mHomopolymerLengths == null)
            return null;

        return mT0AvgQuals.stream()
                .map(x -> format("%.3f", x))
                .collect(Collectors.joining(ITEM_DELIM));
    }
}
