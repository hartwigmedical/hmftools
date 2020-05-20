package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.chaining.SvChainState.EXHAUSTED_PLOIDY_PERC;
import static com.hartwig.hmftools.linx.chaining.SvChainState.EXHAUSTED_PLOIDY_ABS;

public class SegmentPloidy
{
    public final int PosStart;
    public final int PosEnd;
    public final double MajorAP;
    public final double MinorAP;
    public double AFixedAP;
    public double BUndisruptedAP;

    private boolean mIsValid;

    private double mClusterAP;
    private double mLinkPloidy;
    private double mExhaustedLevel;

    public SegmentPloidy(final int posStart, final int posEnd, double major, double minor)
    {
        PosStart = posStart;
        PosEnd = posEnd;
        MajorAP = major;
        MinorAP = minor;
        AFixedAP = 0;
        BUndisruptedAP = 0;
        mClusterAP = 0;

        mIsValid = false;

        mLinkPloidy = 0;
        mExhaustedLevel = 0;
    }

    public boolean isValid() { return mIsValid; }
    public void setValid(boolean toggle) { mIsValid = toggle; }

    public int length() { return PosEnd - PosStart; }
    public double clusterPloidy() { return mClusterAP; }
    public void setClusterPloidy(double ploidy)
    {
        mClusterAP = ploidy;
        mExhaustedLevel = min(EXHAUSTED_PLOIDY_PERC * mClusterAP, EXHAUSTED_PLOIDY_ABS);
    }

    public double unlinkedPloidy() { return max(mClusterAP - mLinkPloidy, 0); }

    public double exhaustedPloidy() { return max(mClusterAP - mLinkPloidy, 0); }

    public void addLinkPloidy(double linkPloidy)
    {
        if(unlinkedPloidy() > linkPloidy)
            mLinkPloidy += linkPloidy;
        else
            mLinkPloidy = mClusterAP;
    }

    public String toString()
    {
        return String.format("%s: A=%.1f B=%.1f C=%.1f Maj=%.1f Min=%.1f Link=%.1f",
                mIsValid ? "valid" : "invalid", AFixedAP, BUndisruptedAP, mClusterAP, MajorAP, MinorAP, mLinkPloidy);
    }
}
