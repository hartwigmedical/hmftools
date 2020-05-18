package com.hartwig.hmftools.common.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstreamTrans;

    private boolean mIsReportable;
    private boolean mPhaseMatched;
    private int mExonsSkippedUp;
    private int mExonsSkippedDown;
    private String mKnownType;
    private boolean mNeoEpitopeOnly;

    private FusionAnnotations mAnnotations;

    // calculated priority accoriding to scheme for selecting fusions
    private double mPriority;

    public static String REPORTABLE_TYPE_NONE = "";
    public static String REPORTABLE_TYPE_KNOWN = "Known";
    public static String REPORTABLE_TYPE_BOTH_PROM = "Both-Prom";
    public static String REPORTABLE_TYPE_5P_PROM = "5P-Prom";
    public static String REPORTABLE_TYPE_3P_PROM = "3P-Prom";

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstream, boolean phaseMatched)
    {
        mUpstreamTrans = upstreamTrans;
        mDownstreamTrans = downstream;
        mIsReportable = false;
        mKnownType = REPORTABLE_TYPE_NONE;
        mPhaseMatched = phaseMatched;
        mExonsSkippedUp = 0;
        mExonsSkippedDown = 0;
        mNeoEpitopeOnly = false;
        mAnnotations = null;
        mPriority = 0;
    }

    public String name() { return mUpstreamTrans.geneName() + "_" + mDownstreamTrans.geneName(); }

    public int svId(boolean isUpstream) { return isUpstream ? mUpstreamTrans.gene().id() : mDownstreamTrans.gene().id(); }

    public Transcript upstreamTrans() { return mUpstreamTrans; }
    public Transcript downstreamTrans() { return mDownstreamTrans; }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public boolean neoEpitopeOnly(){ return mNeoEpitopeOnly; }
    public void setNeoEpitopeOnly(boolean toggle) { mNeoEpitopeOnly = toggle; }

    public final String knownType(){ return mKnownType; }
    public void setKnownType(final String type) { mKnownType = type; }

    public final String toString()
    {
        return String.format("%s %s phased(%s) known(%s)",
                mUpstreamTrans.toString(), mDownstreamTrans.toString(), mPhaseMatched, mKnownType);
    }

    public boolean phaseMatched(){ return mPhaseMatched; }

    public void setExonsSkipped(int exonsUp, int exonsDown)
    {
        mExonsSkippedUp = exonsUp;
        mExonsSkippedDown = exonsDown;
    }

    public int getExonsSkipped(boolean isUpstream) { return isUpstream ? mExonsSkippedUp : mExonsSkippedDown; }

    public int getFusedExon(boolean isUpstream)
    {
        if(isUpstream)
        {
            return max(mUpstreamTrans.ExonUpstream - mExonsSkippedUp, 1);
        }
        else
        {
            return min(mDownstreamTrans.ExonDownstream + mExonsSkippedDown, mDownstreamTrans.ExonMax);
        }
    }

    public boolean isExonic()
    {
        return mUpstreamTrans.isExonic() &&mDownstreamTrans.isExonic();
    }

    public final FusionAnnotations getAnnotations() { return mAnnotations; }
    public void setAnnotations(final FusionAnnotations annotations) { mAnnotations = annotations; }

    public void setPriority(double priority) { mPriority = priority; }
    public double priority() { return mPriority; }

    public boolean isTerminated()
    {
        if(mAnnotations == null)
            return false;

        if((mAnnotations.disruptionUp() != null && mAnnotations.disruptionUp().transcriptTerminated())
        || (mAnnotations.disruptionDown() != null && mAnnotations.disruptionDown().transcriptTerminated()))
        {
            return true;
        }

        return false;
    }

    // convenience functions
    public int getChainLength()
    {
        return mAnnotations != null && mAnnotations.chainInfo() != null ? mAnnotations.chainInfo().chainLength() : 0;
    }

    public int getChainLinks()
    {
        return mAnnotations != null && mAnnotations.chainInfo() != null ? mAnnotations.chainInfo().chainLinks() : 0;
    }

    public boolean validChainTraversal()
    {
        if(mAnnotations != null && mAnnotations.chainInfo() != null)
            return mAnnotations.chainInfo().validTraversal();
        else
            return true;
    }

    public boolean isViable()
    {
        if(!validChainTraversal() || isTerminated())
            return false;

        if(mDownstreamTrans.hasNegativePrevSpliceAcceptorDistance())
            return false;

        return true;
    }

}
