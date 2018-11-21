package com.hartwig.hmftools.common.variant.structural.annotation;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstream;

    private final String mPrimarySource;

    private final boolean mIsPhaseMatch;

    private boolean mIsReportable;

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstream, final String primarySource, boolean isReportable, boolean isPhaseMatch)
    {
        mUpstreamTrans = upstreamTrans;
        mDownstream = downstream;
        mPrimarySource = primarySource;
        mIsReportable = isReportable;
        mIsPhaseMatch = isPhaseMatch;
    }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public Transcript upstreamTrans() { return mUpstreamTrans; }

    public Transcript downstreamTrans() { return mDownstream; }

    public String primarySource() { return mPrimarySource; }

    public boolean isPhaseMatch() { return mIsPhaseMatch; }
}
