package com.hartwig.hmftools.common.basequal.jitter;

public class JitterModelLoss
{
    private JitterCountsTable.Row mStatsTableRow;
    private double mNumRepeats;
    private double mLengthMinusOneScale;

    public JitterModelLoss(JitterCountsTable.Row statsTableRow, double numRepeats, double lengthMinusOneScale)
    {
        this.mStatsTableRow = statsTableRow;
        this.mNumRepeats = numRepeats;
        this.mLengthMinusOneScale = lengthMinusOneScale;
    }

    // total loss is the sum of loss from jitter values from -5 to 5
    public double totalLoss(double scale, double skew)
    {
        double lossSum = 0;
        for(int x = -5; x <= 5; ++x)
        {
            if(x != 0)
            {
                lossSum += loss(x, scale, skew);
            }
        }

        lossSum += (Math.max(skew, 1.0 / skew) - 1) * 0.01;

        if(mNumRepeats >= 7)
        {
            lossSum += Math.pow(Math.abs(scale - mLengthMinusOneScale), 2);
        }

        return lossSum;
    }

    public double loss(int x, double scale, double skew)
    {
        double modelFreq = Math.max(JitterModelCalc.modifiedAsymmetricLaplace(x, scale, skew), 1e-4);
        double flooredFreq = Math.max(freq(x), 1e-4);
        return -10 * Math.log10(Math.min(flooredFreq, modelFreq) / Math.max(flooredFreq, modelFreq)) * freq(x);
    }

    private double freq(int x)
    {
        return ((double)mStatsTableRow.getJitterReadCount(x)) / mStatsTableRow.totalReadCount;
    }
}
