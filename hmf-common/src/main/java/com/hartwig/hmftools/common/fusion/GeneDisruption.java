package com.hartwig.hmftools.common.fusion;

public class GeneDisruption
{
    private final Transcript mTranscript;
    private boolean mIsReportable;

    public GeneDisruption(final Transcript transcript)
    {
        mTranscript = transcript;
        mIsReportable = false;
    }

    public final Transcript transcript() { return mTranscript; }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

}
