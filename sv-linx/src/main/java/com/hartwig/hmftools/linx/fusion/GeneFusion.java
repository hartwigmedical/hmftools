package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.UPSTREAM_STR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.fsIndex;

import com.hartwig.hmftools.common.fusion.Transcript;

public class GeneFusion
{
    private final Transcript[] mTranscripts;
    private final String mName;

    private boolean mIsReportable;
    private boolean mPhaseMatched;
    private int[] mExonsSkipped;
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

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstreamTrans, boolean phaseMatched)
    {
        mTranscripts = new Transcript[] { upstreamTrans, downstreamTrans };
        mName = mTranscripts[FS_UPSTREAM].geneName() + "_" + mTranscripts[FS_DOWNSTREAM].geneName();
        mIsReportable = false;
        mKnownType = REPORTABLE_TYPE_NONE;
        mPhaseMatched = phaseMatched;
        mExonsSkipped = new int[] { 0, 0 };
        mNeoEpitopeOnly = false;
        mAnnotations = null;
        mPriority = 0;
    }

    public String name() { return mName; }

    public int svId(boolean isUpstream) { return mTranscripts[fsIndex(isUpstream)].gene().id(); }

    public Transcript transcript(int fs) { return mTranscripts[fs]; }
    public Transcript upstreamTrans() { return mTranscripts[FS_UPSTREAM]; }
    public Transcript downstreamTrans() { return mTranscripts[FS_DOWNSTREAM]; }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public boolean neoEpitopeOnly(){ return mNeoEpitopeOnly; }
    public void setNeoEpitopeOnly(boolean toggle) { mNeoEpitopeOnly = toggle; }

    public final String knownType(){ return mKnownType; }
    public void setKnownType(final String type) { mKnownType = type; }

    public final String toString()
    {
        return String.format("%s %s phased(%s) known(%s)",
                mTranscripts[FS_UPSTREAM].toString(), mTranscripts[FS_DOWNSTREAM].toString(), mPhaseMatched, mKnownType);
    }

    public boolean phaseMatched(){ return mPhaseMatched; }

    public void setExonsSkipped(int exonsUp, int exonsDown)
    {
        mExonsSkipped[FS_UPSTREAM] = exonsUp;
        mExonsSkipped[FS_DOWNSTREAM] = exonsDown;
    }

    public int getExonsSkipped(boolean isUpstream) { return mExonsSkipped[fsIndex(isUpstream)]; }
    public int[] getExonsSkipped() { return mExonsSkipped; }

    public int getFusedExon(boolean isUpstream)
    {
        if(isUpstream)
        {
            return max(mTranscripts[FS_UPSTREAM].ExonUpstream - mExonsSkipped[FS_UPSTREAM], 1);
        }
        else
        {
            return min(mTranscripts[FS_DOWNSTREAM].ExonDownstream + mExonsSkipped[FS_DOWNSTREAM], mTranscripts[FS_DOWNSTREAM].ExonMax);
        }
    }

    public boolean isExonic()
    {
        return mTranscripts[FS_UPSTREAM].isExonic() && mTranscripts[FS_DOWNSTREAM].isExonic();
    }

    public final FusionAnnotations getAnnotations() { return mAnnotations; }
    public void setAnnotations(final FusionAnnotations annotations) { mAnnotations = annotations; }

    public void setPriority(double priority) { mPriority = priority; }
    public double priority() { return mPriority; }

    public boolean isTerminated()
    {
        return mAnnotations != null && (mAnnotations.terminatedUp() || mAnnotations.terminatedDown());
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

        if(mTranscripts[FS_DOWNSTREAM].hasNegativePrevSpliceAcceptorDistance())
            return false;

        return true;
    }

}
