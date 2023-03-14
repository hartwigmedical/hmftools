package com.hartwig.hmftools.ctdna.purity;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class CopyNumberGcData
{
    public final String Chromosome;
    public final int SegmentStart;
    public final int SegmentEnd;
    public final double CopyNumber;

    private final List<Double> mTumorGcRatios;

    private boolean mComputed;
    private double mMean;
    private double mMedian;

    public CopyNumberGcData(
            final String chromosome, final int segmentStart, final int segmentEnd, final double copyNumber)
    {
        Chromosome = chromosome;
        SegmentStart = segmentStart;
        SegmentEnd = segmentEnd;
        CopyNumber = copyNumber;

        mTumorGcRatios = Lists.newArrayList();

        mComputed = false;
        mMean = 0;
        mMedian = 0;
    }

    public void addRatio(double ratio)
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

        Collections.sort(mTumorGcRatios);

        int index = mTumorGcRatios.size() / 2;

        if((mTumorGcRatios.size() % 2) == 0)
        {
            mMedian = (mTumorGcRatios.get(index - 1) + mTumorGcRatios.get(index)) * 0.5;
        }
        else
        {
            mMedian = mTumorGcRatios.get(index);
        }

        double total = mTumorGcRatios.stream().mapToDouble(x -> x).sum();
        mMean = total / mTumorGcRatios.size();

        mComputed = true;
    }
}
