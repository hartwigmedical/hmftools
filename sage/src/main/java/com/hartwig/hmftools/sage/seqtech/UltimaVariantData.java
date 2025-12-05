package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.HALF_PHRED_SCORE_SCALING;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractTpValues;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;
import static com.hartwig.hmftools.sage.vcf.ReadContextVcfInfo.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class UltimaVariantData
{
    private final VariantReadContext mReadContext;

    private UltimaRealignedQualModels mQualModels;

    private final List<Integer> mHomopolymerLengths;

    private final int mHomopolymerLeftPadding;
    private final int mHomopolymerRightPadding;

    private final List<Double> mHomopolymerAvgQuals;
    private final List<Double> mT0AvgQuals;
    private int mHomopolymerAvgQualsCount;

    public UltimaVariantData(final VariantReadContext readContext)
    {
        mReadContext = readContext;

        List<Homopolymer> homopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);

        mHomopolymerLengths = homopolymers.stream().map(x -> x.Length).collect(Collectors.toList());

        mHomopolymerLeftPadding = getHomopolymerPadding(homopolymers, readContext, true);
        mHomopolymerRightPadding = getHomopolymerPadding(homopolymers, readContext, false);

        mHomopolymerAvgQuals = Lists.newArrayList();
        mT0AvgQuals = Lists.newArrayList();

        while(mHomopolymerAvgQuals.size() < mHomopolymerLengths.size())
        {
            mHomopolymerAvgQuals.add(0.0);
            mT0AvgQuals.add(0.0);
        }
    }

    public List<Integer> paddedHomopolymerLengths()
    {
        // add HP extensions from the flanks to the outer homopolymers
        List<Integer> paddedLengths = Lists.newArrayList();

        for(int i = 0; i < mHomopolymerLengths.size(); ++i)
        {
            int hpLength = mHomopolymerLengths.get(i);

            if(i == 0)
                hpLength += mHomopolymerLeftPadding;

            if(i == mHomopolymerLengths.size() - 1)
                hpLength += mHomopolymerRightPadding;

            paddedLengths.add(hpLength);
        }

        return paddedLengths;
    }

    public List<Double> homopolymerAvgQuals() { return mHomopolymerAvgQuals; }
    public List<Double> t0AvgQuals() { return mT0AvgQuals; }

    public UltimaRealignedQualModels getQualModels() { return mQualModels; }
    public void setQualModels(final UltimaRealignedQualModels models) { mQualModels = models; }

    public void addReadSupportInfo(final SAMRecord record, final int readVarIndex)
    {
        // do an exact match since quals are not applicable
        String recordCore = record.getReadString().substring(
                readVarIndex - mReadContext.leftCoreLength(), readVarIndex + mReadContext.rightCoreLength() + 1);

        if(recordCore.equals(mReadContext.coreStr()))
        {
            registerHomopolymerQuals(record, readVarIndex);
            registerT0Quals(record, readVarIndex);
        }
    }

    private int getHomopolymerPadding(final List<Homopolymer> homopolymers, final VariantReadContext readContext, boolean isLeft)
    {
        int hpPadding = 0;

        if(isLeft)
        {
            for(int i = readContext.CoreIndexStart - 1; i >= 0; --i)
            {
                if(readContext.ReadBases[i] == homopolymers.get(0).Base)
                    ++hpPadding;
                else
                    break;
            }
        }
        else
        {
            for(int i = readContext.CoreIndexEnd + 1; i < readContext.ReadBases.length; ++i)
            {
                if(readContext.ReadBases[i] == homopolymers.get(homopolymers.size() - 1).Base)
                    ++hpPadding;
                else
                    break;
            }
        }

        return hpPadding;
    }

    private void registerHomopolymerQuals(final SAMRecord record, int readVarIndex)
    {
        byte[] baseQuals = record.getBaseQualities();
        byte[] tpValues = extractTpValues(record);
        int readIndex = readVarIndex - mReadContext.leftCoreLength();
        boolean firstHomopolymer = true;
        List<Byte> homopolyerQuals = Lists.newArrayList();

        for(int hpLength : mHomopolymerLengths)
        {
            int lookupIndex = firstHomopolymer ? readIndex + hpLength - 1 : readIndex;
            int tpValue = tpValues[lookupIndex];

            byte homopolymerQual = baseQuals[lookupIndex];

            if(tpValue != 0 && hpLength > 1)
                homopolymerQual -= HALF_PHRED_SCORE_SCALING;

            homopolyerQuals.add(homopolymerQual);

            firstHomopolymer = false;
            readIndex += hpLength;
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

        List<Byte> t0Quals = Lists.newArrayList();

        for(int hpLength : mHomopolymerLengths)
        {
            int lookupIndex = firstHomopolymer ? readIndex + hpLength - 1 : readIndex;
            byte t0Value = t0Values[lookupIndex];

            t0Quals.add(t0Value);

            firstHomopolymer = false;
            readIndex += hpLength;
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

        String lengthsStr = paddedHomopolymerLengths().stream()
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
