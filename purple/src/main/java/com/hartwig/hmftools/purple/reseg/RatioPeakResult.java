package com.hartwig.hmftools.purple.reseg;

public final class RatioPeakResult
{
    public final PeakTroughData PrimaryPeak;
    public final double LeftBound;
    public final double RightBound;

    public RatioPeakResult(final PeakTroughData primaryPeak, final double leftBound, final double rightBound)
    {
        PrimaryPeak = primaryPeak;
        LeftBound = leftBound;
        RightBound = rightBound;
    }

    public String toString()
    {
        return String.format("primaryPeak(%s) bounds(%.4f - %.4f)", PrimaryPeak, LeftBound, RightBound);
    }
}
