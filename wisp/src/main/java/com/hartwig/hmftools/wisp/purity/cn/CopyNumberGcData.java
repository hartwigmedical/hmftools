package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class CopyNumberGcData
{
    public final String Chromosome;
    public final int SegmentStart;
    public final int SegmentEnd;
    public final double CopyNumber;

    public final int CopyNumberLevel;
    public final boolean IsValid;

    private final List<GcRatioData> mTumorGcRatios;

    private boolean mComputed;
    private double mMean;
    private double mMedian;

    public CopyNumberGcData(
            final String chromosome, final int segmentStart, final int segmentEnd, final double copyNumber,
            boolean isValid, int copyNumberLevel)
    {
        Chromosome = chromosome;
        SegmentStart = segmentStart;
        SegmentEnd = segmentEnd;
        CopyNumber = copyNumber;
        CopyNumberLevel = copyNumberLevel;
        IsValid = isValid;

        mTumorGcRatios = Lists.newArrayList();

        mComputed = false;
        mMean = 0;
        mMedian = 0;
    }

    public List<GcRatioData> ratios() { return mTumorGcRatios; }

    public void addRatio(final GcRatioData ratio)
    {
        mTumorGcRatios.add(ratio);
        mComputed = false;
    }

    public int count() { return mTumorGcRatios.size(); }

    public double mean()
    {
        if(!mComputed)
            compute();

        return mMean;
    }

    public double median()
    {
        if(!mComputed)
            compute();

        return mMedian;
    }

    private void compute()
    {
        if(mTumorGcRatios.isEmpty())
        {
            mMedian = 0;
            mMean = 0;
            mComputed = true;
            return;
        }

        List<GcRatioData> sortedList = Lists.newArrayList(mTumorGcRatios);
        Collections.sort(sortedList);

        int index = sortedList.size() / 2;

        if((sortedList.size() % 2) == 0)
        {
            mMedian = (sortedList.get(index - 1).TumorGcRatio + sortedList.get(index).TumorGcRatio) * 0.5;
        }
        else
        {
            mMedian = sortedList.get(index).TumorGcRatio;
        }

        double total = sortedList.stream().mapToDouble(x -> x.TumorGcRatio).sum();
        mMean = total / sortedList.size();

        mComputed = true;
    }

    public String toString()
    {
        String coreData = format("loc(%s:%d-%d) cn(%.2f) ratios(%d)", Chromosome, SegmentStart, SegmentEnd, CopyNumber, mTumorGcRatios.size());

        if(!mComputed)
            return coreData;

        return format("%s mean(%.4f) median(%.4f)", coreData, mMean, mMedian);
    }
}
