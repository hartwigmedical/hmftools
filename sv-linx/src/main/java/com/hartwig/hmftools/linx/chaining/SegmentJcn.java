package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.chaining.ChainState.EXHAUSTED_JCN_PERC;
import static com.hartwig.hmftools.linx.chaining.ChainState.EXHAUSTED_JCN_ABS;

public class SegmentJcn
{
    public final int PosStart;
    public final int PosEnd;
    public final double MajorAP;
    public final double MinorAP;
    public double AFixedAP;
    public double BUndisruptedAP;

    private boolean mIsValid;

    private double mClusterAP;
    private double mLinkJcn;
    private double mExhaustedLevel;

    public SegmentJcn(final int posStart, final int posEnd, double major, double minor)
    {
        PosStart = posStart;
        PosEnd = posEnd;
        MajorAP = major;
        MinorAP = minor;
        AFixedAP = 0;
        BUndisruptedAP = 0;
        mClusterAP = 0;

        mIsValid = false;

        mLinkJcn = 0;
        mExhaustedLevel = 0;
    }

    public boolean isValid() { return mIsValid; }
    public void setValid(boolean toggle) { mIsValid = toggle; }

    public int length() { return PosEnd - PosStart; }
    public double clusterJcn() { return mClusterAP; }
    public void setClusterJcn(double jcn)
    {
        mClusterAP = jcn;
        mExhaustedLevel = min(EXHAUSTED_JCN_PERC * mClusterAP, EXHAUSTED_JCN_ABS);
    }

    public double unlinkedJcn() { return max(mClusterAP - mLinkJcn, 0); }

    public double exhaustedJcn() { return max(mClusterAP - mLinkJcn, 0); }

    public void addLinkJcn(double linkPloidy)
    {
        if(unlinkedJcn() > linkPloidy)
            mLinkJcn += linkPloidy;
        else
            mLinkJcn = mClusterAP;
    }

    public String toString()
    {
        return String.format("%s: A=%.1f B=%.1f C=%.1f Maj=%.1f Min=%.1f Link=%.1f",
                mIsValid ? "valid" : "invalid", AFixedAP, BUndisruptedAP, mClusterAP, MajorAP, MinorAP, mLinkJcn);
    }
}
