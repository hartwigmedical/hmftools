package com.hartwig.hmftools.common.variant.structural.annotation;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstreamTrans;

    private boolean mIsReportable;
    private boolean mPhaseMatched;
    private boolean mViable; // passes fusion rules
    private int mExonsSkippedUp;
    private int mExonsSkippedDown;
    private String mKnownType;
    private boolean mNeoEpitopeOnly;

    // private String mAnnotations;
    private FusionAnnotations mAnnotations;

    public static String REPORTABLE_TYPE_NONE = "";
    public static String REPORTABLE_TYPE_KNOWN = "Known";
    public static String REPORTABLE_TYPE_BOTH_PROM = "Both-Prom";
    public static String REPORTABLE_TYPE_5P_PROM = "5P-Prom";
    public static String REPORTABLE_TYPE_3P_PROM = "3P-Prom";

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstream, boolean phaseMatched, boolean viable)
    {
        mUpstreamTrans = upstreamTrans;
        mDownstreamTrans = downstream;
        mIsReportable = false;
        mKnownType = REPORTABLE_TYPE_NONE;
        mPhaseMatched = phaseMatched;
        mViable = viable;
        mExonsSkippedUp = 0;
        mExonsSkippedDown = 0;
        mNeoEpitopeOnly = false;
        mAnnotations = null;
    }

    public String name() { return mUpstreamTrans.geneName() + "_" + mDownstreamTrans.geneName(); }

    public int svId(boolean isUpstream) { return isUpstream ? mUpstreamTrans.gene().id() : mDownstreamTrans.gene().id(); }

    public Transcript upstreamTrans() { return mUpstreamTrans; }
    public Transcript downstreamTrans() { return mDownstreamTrans; }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public boolean neoEpitopeOnly(){ return mNeoEpitopeOnly; }
    public void setNeoEpitopeOnly(boolean toggle) { mNeoEpitopeOnly = toggle; }

    public final String getKnownType(){ return mKnownType; }
    public void setKnownType(final String type) { mKnownType = type; }

    public boolean phaseMatched(){ return mPhaseMatched; }
    public boolean viable(){ return mViable; }
    public void setViable(boolean viable){ mViable = viable; }

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
            if (mUpstreamTrans.isExonic() && !mDownstreamTrans.isExonic())
                return mUpstreamTrans.ExonUpstream - 1 - mExonsSkippedUp;
            else
                return mUpstreamTrans.ExonUpstream - mExonsSkippedUp;
        }
        else
        {
            if (mDownstreamTrans.isExonic() && !mUpstreamTrans.isExonic())
                return mDownstreamTrans.ExonDownstream + 1 + mExonsSkippedDown;
            else
                return mDownstreamTrans.ExonDownstream + mExonsSkippedDown;
        }
    }

    public boolean isExonic()
    {
        return mUpstreamTrans.isExonic() &&mDownstreamTrans.isExonic();
    }

    public final FusionAnnotations getAnnotations() { return mAnnotations; }
    public void setAnnotations(final FusionAnnotations annotations) { mAnnotations = annotations; }

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
    public long getChainLength()
    {
        return mAnnotations != null && mAnnotations.chainInfo() != null ? mAnnotations.chainInfo().chainLength() : 0;
    }

    public long getChainLinks()
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
        return mPhaseMatched && validChainTraversal() && !isTerminated() && mDownstreamTrans.exonDistanceUp() > 0;
    }

}
