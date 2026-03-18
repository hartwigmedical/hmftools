package com.hartwig.hmftools.redux.ms_model;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.Map;

import com.hartwig.hmftools.common.redux.JitterTableRow;

public class RepeatUnitData
{
    public final String RepeatUnit;
    public final JitterTableRow JitterRow;

    private int mWeightedErrorCount;
    private int mWeightedTotal;
    private double mErrorRate;
    private double mAdjustedErrorRate;

    private final String mRucKey;

    public RepeatUnitData(final String repeatUnit, final JitterTableRow jitterRow)
    {
        RepeatUnit = repeatUnit;
        JitterRow = jitterRow;
        mRucKey = repeatUnitAndCountKey(RepeatUnit, jitterRow.refNumUnits());

        mAdjustedErrorRate = 0;
    }

    public String rucKey()
    {
        return mRucKey;
    }

    public int repeatCount()
    {
        return JitterRow.refNumUnits();
    }

    public void mergeJitterRow(final JitterTableRow jcRow)
    {
        for(Map.Entry<Integer, Integer> entry : jcRow.jitterCounts().entrySet())
        {
            JitterRow.addReads(entry.getKey(), entry.getValue());
        }
    }

    public void computeErrorRates()
    {
        mWeightedErrorCount = 0;
        mWeightedTotal = 0;
        mErrorRate = 0;

        for(Map.Entry<Integer, Integer> entry : JitterRow.jitterCounts().entrySet())
        {
            int count = abs(entry.getKey());
            int reads = entry.getValue();

            if(count == 0)
            {
                mWeightedTotal += reads;
            }
            else
            {
                mWeightedErrorCount += count * reads;
            }
        }

        mWeightedTotal += mWeightedErrorCount;
        mErrorRate = mWeightedTotal > 0 ? mWeightedErrorCount / (double) mWeightedTotal : 0;
    }

    public void setAdjustedErrorRate(double noiseAdjustment)
    {
        mAdjustedErrorRate = max(mErrorRate - noiseAdjustment, 0);
    }
    public double adjustedErrorRate()
    {
        return mAdjustedErrorRate;
    }

    public double errorRate()
    {
        return mErrorRate;
    }

    public String toString()
    {
        return format("%s: reads(%d) errorRate(%4.3e)", mRucKey, mWeightedTotal, mErrorRate);
    }

    protected static String repeatUnitAndCountKey(final String repeatUnit, final int repeatCount)
    {
        // eg A/T with num-units = 4
        return format("%s_%d", repeatUnit, repeatCount);
    }

}
