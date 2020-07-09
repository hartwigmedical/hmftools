package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.UPSTREAM_STR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.fsIndex;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_BOTH;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.Transcript;

public class GeneFusion
{
    private final Transcript[] mTranscripts;

    private boolean mIsReportable;
    private boolean mPhaseMatched;
    private int[] mExonsSkipped;
    private KnownFusionType mKnownFusionType;
    private final boolean[] mIsPromiscuous;
    private boolean mNeoEpitopeOnly;

    private FusionAnnotations mAnnotations;

    // calculated priority accoriding to scheme for selecting fusions
    private double mPriority;

    public GeneFusion(final Transcript upstreamTrans, final Transcript downstreamTrans, boolean phaseMatched)
    {
        mTranscripts = new Transcript[] { upstreamTrans, downstreamTrans };

        mIsReportable = false;
        mKnownFusionType = KnownFusionType.NONE;
        mIsPromiscuous = new boolean[] { false, false };
        mPhaseMatched = phaseMatched;
        mExonsSkipped = new int[] { 0, 0 };
        mNeoEpitopeOnly = false;
        mAnnotations = null;
        mPriority = 0;
    }

    public String name()
    {
        if(mKnownFusionType == IG_KNOWN_PAIR || mKnownFusionType == IG_PROMISCUOUS) // form like @IGH-MYC
            return mTranscripts[FS_UPSTREAM].StableId + "_" + mTranscripts[FS_DOWNSTREAM].geneName();
        else
            return mTranscripts[FS_UPSTREAM].geneName() + "_" + mTranscripts[FS_DOWNSTREAM].geneName();
    }

    public int svId(boolean isUpstream) { return mTranscripts[fsIndex(isUpstream)].gene().id(); }

    public Transcript[] transcripts() { return mTranscripts; }
    public Transcript upstreamTrans() { return mTranscripts[FS_UPSTREAM]; }
    public Transcript downstreamTrans() { return mTranscripts[FS_DOWNSTREAM]; }

    public String geneName(int fs)
    {
        if(fs == FS_UPSTREAM && (mKnownFusionType == IG_KNOWN_PAIR || mKnownFusionType == IG_PROMISCUOUS))
            return mTranscripts[FS_UPSTREAM].StableId;
        else
            return mTranscripts[fs].geneName();
    }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }

    public boolean neoEpitopeOnly(){ return mNeoEpitopeOnly; }
    public void setNeoEpitopeOnly(boolean toggle) { mNeoEpitopeOnly = toggle; }

    public final String knownTypeStr()
    {
        if(mKnownFusionType == PROMISCUOUS_5 || mKnownFusionType == PROMISCUOUS_3)
        {
            if(mIsPromiscuous[FS_UPSTREAM] && mIsPromiscuous[FS_DOWNSTREAM])
                return PROMISCUOUS_BOTH;
        }

        return mKnownFusionType.toString();
    }

    public KnownFusionType knownType() { return mKnownFusionType; }
    public boolean[] isPromiscuous() { return mIsPromiscuous; }

    public void setKnownType(KnownFusionType type)
    {
        mKnownFusionType = type;
    }

    public boolean phaseMatched(){ return mPhaseMatched; }

    public void setExonsSkipped(int exonsUp, int exonsDown)
    {
        mExonsSkipped[FS_UPSTREAM] = exonsUp;
        mExonsSkipped[FS_DOWNSTREAM] = exonsDown;
    }

    public boolean sameChromosome()
    {
        return mTranscripts[FS_UPSTREAM].gene().chromosome().equals(mTranscripts[FS_DOWNSTREAM].gene().chromosome());
    }

    public int distance()
    {
        if(!sameChromosome())
            return 0;

        return abs(mTranscripts[FS_UPSTREAM].gene().position() - mTranscripts[FS_DOWNSTREAM].gene().position());
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
        if(!validChainTraversal() || (isTerminated() && mKnownFusionType != KNOWN_PAIR))
            return false;

        if(mTranscripts[FS_DOWNSTREAM].hasNegativePrevSpliceAcceptorDistance())
            return false;

        return true;
    }

    public final String toString()
    {
        return String.format("%s %s phased(%s) type(%s)",
                mTranscripts[FS_UPSTREAM].toString(), mTranscripts[FS_DOWNSTREAM].toString(), mPhaseMatched, knownTypeStr());
    }
}
