package com.hartwig.hmftools.common.variant.structural.annotation;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstream;

    private final String mPrimarySource;

    private final boolean mIsPhaseMatch;

    private boolean mIsReportable;

    public static String RNA_MATCH_UNKNOWN = "Unknown";
    public static String RNA_MATCH_MATCHED = "Match";
    public static String RNA_MATCH_NO_RNA = "NoRNA";
    public static String RNA_MATCH_NO_SV_FUSION = "NoSVFusion";
    private String mRnaMatch;

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstream, final String primarySource, boolean isReportable, boolean isPhaseMatch)
    {
        mUpstreamTrans = upstreamTrans;
        mDownstream = downstream;
        mPrimarySource = primarySource;
        mIsReportable = isReportable;
        mIsPhaseMatch = isPhaseMatch;
        mRnaMatch = RNA_MATCH_UNKNOWN;
    }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public final String getRnaMatchType(){ return mRnaMatch; }
    public void setRnaMatch(final String type) { mRnaMatch = type; }

    public Transcript upstreamTrans() { return mUpstreamTrans; }

    public Transcript downstreamTrans() { return mDownstream; }

    public String primarySource() { return mPrimarySource; }

    public boolean isPhaseMatch() { return mIsPhaseMatch; }
}
