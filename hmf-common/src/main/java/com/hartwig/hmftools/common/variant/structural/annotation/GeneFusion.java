package com.hartwig.hmftools.common.variant.structural.annotation;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstreamTrans;

    private String mPrimarySource;

    private boolean mIsReportable;
    private boolean mPhaseMatched;
    private boolean mViable; // passes fusion rules
    private int mExonsSkippedUp;
    private int mExonsSkippedDown;
    private String mKnownFusionType;

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
        mPrimarySource = "";
        mIsReportable = false;
        mKnownFusionType = REPORTABLE_TYPE_NONE;
        mPhaseMatched = phaseMatched;
        mViable = viable;
        mExonsSkippedUp = 0;
        mExonsSkippedDown = 0;
        mAnnotations = null;
    }

    public String name() { return mUpstreamTrans.geneName() + "_" + mDownstreamTrans.geneName(); }

    public Transcript upstreamTrans() { return mUpstreamTrans; }
    public Transcript downstreamTrans() { return mDownstreamTrans; }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public final String getKnownFusionType(){ return mKnownFusionType; }
    public void setKnownFusionType(final String type) { mKnownFusionType = type; }

    public String primarySource() { return mPrimarySource; }
    public void setPrimarySource(final String source) { mPrimarySource = source; }

    public boolean phaseMatched(){ return mPhaseMatched; }
    public boolean viable(){ return mViable; }
    public void setExonsSkipped(int exonsUp, int exonsDown)
    {
        mExonsSkippedUp = exonsUp;
        mExonsSkippedDown = exonsDown;
    }

    public int getExonsSkipped(boolean isUpstream) { return isUpstream ? mExonsSkippedUp : mExonsSkippedDown; }

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
        return mPhaseMatched && mExonsSkippedUp == 0 && mExonsSkippedDown == 0 && validChainTraversal() && !isTerminated();
    }

}
