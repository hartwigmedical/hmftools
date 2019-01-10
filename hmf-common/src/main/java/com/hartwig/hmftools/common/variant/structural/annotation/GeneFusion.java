package com.hartwig.hmftools.common.variant.structural.annotation;

public class GeneFusion
{
    private final Transcript mUpstreamTrans;

    private final Transcript mDownstream;

    private final String mPrimarySource;

    private boolean mIsReportable;
    private String mKnownFusionType;

    public static String REPORTABLE_TYPE_NONE = "";
    public static String REPORTABLE_TYPE_KNOWN = "Known";
    public static String REPORTABLE_TYPE_BOTH_PROM = "Both-Prom";
    public static String REPORTABLE_TYPE_5P_PROM = "5P-Prom";
    public static String REPORTABLE_TYPE_3P_PROM = "3P-Prom";

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstream, final String primarySource, boolean isReportable)
    {
        mUpstreamTrans = upstreamTrans;
        mDownstream = downstream;
        mPrimarySource = primarySource;
        mIsReportable = isReportable;
        mKnownFusionType = REPORTABLE_TYPE_NONE;
    }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public final String getKnownFusionType(){ return mKnownFusionType; }
    public void setKnownFusionType(final String type) { mKnownFusionType = type; }

    public Transcript upstreamTrans() { return mUpstreamTrans; }

    public Transcript downstreamTrans() { return mDownstream; }

    public String primarySource() { return mPrimarySource; }
}
