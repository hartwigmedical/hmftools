package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.fsIndex;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_BOTH;
import static com.hartwig.hmftools.common.sv.linx.FusionPhasedType.INFRAME;
import static com.hartwig.hmftools.common.sv.linx.FusionPhasedType.OUT_OF_FRAME;
import static com.hartwig.hmftools.common.sv.linx.FusionPhasedType.SKIPPED_EXONS;
import static com.hartwig.hmftools.linx.fusion.ReportableReason.UNSET;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;

public class GeneFusion
{
    private int mId; // optional identifier
    private final BreakendTransData[] mTranscripts;

    private boolean mIsReportable;
    private ReportableReason mReportableReason;
    private boolean mPhaseMatched;
    private int[] mExonsSkipped;
    private KnownFusionType mKnownFusionType;
    private final boolean[] mIsPromiscuous;
    private boolean mKnownExons;
    private boolean mHighImpactPromiscuous;
    private boolean mProteinFeaturesSet;

    private FusionAnnotations mAnnotations;

    // calculated priority according to scheme for selecting fusions
    private double mPriority;

    public GeneFusion(final BreakendTransData upstreamTrans, final BreakendTransData downstreamTrans, boolean phaseMatched)
    {
        mTranscripts = new BreakendTransData[] { upstreamTrans, downstreamTrans };

        mIsReportable = false;
        mReportableReason = UNSET;
        mKnownFusionType = KnownFusionType.NONE;
        mIsPromiscuous = new boolean[] { false, false };
        mPhaseMatched = phaseMatched;
        mExonsSkipped = new int[] { 0, 0 };
        mKnownExons = false;
        mHighImpactPromiscuous = false;
        mProteinFeaturesSet = false;
        mAnnotations = null;
        mPriority = 0;
    }

    public void setId(final int id) { mId = id; }
    public int id() { return mId; }

    public String name()
    {
        if(mKnownFusionType == IG_KNOWN_PAIR || mKnownFusionType == IG_PROMISCUOUS) // form like @IGH-MYC
            return mTranscripts[FS_UP].transName() + "_" + mTranscripts[FS_DOWN].geneName();
        else
            return mTranscripts[FS_UP].geneName() + "_" + mTranscripts[FS_DOWN].geneName();
    }

    public int svId(boolean isUpstream) { return mTranscripts[fsIndex(isUpstream)].gene().id(); }

    public BreakendTransData[] transcripts() { return mTranscripts; }
    public BreakendTransData upstreamTrans() { return mTranscripts[FS_UP]; }
    public BreakendTransData downstreamTrans() { return mTranscripts[FS_DOWN]; }

    public String geneName(int fs)
    {
        if(fs == FS_UP && (mKnownFusionType == IG_KNOWN_PAIR || mKnownFusionType == IG_PROMISCUOUS))
            return mTranscripts[FS_UP].transName();
        else
            return mTranscripts[fs].geneName();
    }

    public boolean reportable(){ return mIsReportable; }
    public void setReportable(boolean toggle) { mIsReportable = toggle; }
    public void setReportableReason(final ReportableReason reason) { mReportableReason = reason; }
    public ReportableReason reportableReason() { return mReportableReason; }

    public boolean proteinFeaturesSet() { return mProteinFeaturesSet; }
    public void setProteinFeaturesSet() { mProteinFeaturesSet = true; }

    public final String knownTypeStr()
    {
        if(mKnownFusionType == PROMISCUOUS_5 || mKnownFusionType == PROMISCUOUS_3)
        {
            if(mIsPromiscuous[FS_UP] && mIsPromiscuous[FS_DOWN])
                return PROMISCUOUS_BOTH;
        }

        return mKnownFusionType.toString();
    }

    public KnownFusionType knownType() { return mKnownFusionType; }
    public boolean[] isPromiscuous() { return mIsPromiscuous; }

    public void setKnownType(KnownFusionType type) { mKnownFusionType = type; }

    public boolean phaseMatched() { return mPhaseMatched; }

    public FusionPhasedType phaseType()
    {
        if(!mPhaseMatched)
            return OUT_OF_FRAME;

        return mExonsSkipped[FS_UP] > 0 || mExonsSkipped[FS_DOWN] > 0 ? SKIPPED_EXONS : INFRAME;
    }

    public FusionLikelihoodType likelihoodType()
    {
        if(!mIsReportable)
            return FusionLikelihoodType.NA;

        if(mKnownFusionType == KNOWN_PAIR || mKnownFusionType == EXON_DEL_DUP || mKnownFusionType == IG_KNOWN_PAIR || mKnownExons)
            return FusionLikelihoodType.HIGH;
        else
            return FusionLikelihoodType.LOW;
    }

    public void setKnownExons() { mKnownExons = true; }
    public boolean knownExons() { return mKnownExons; }

    public void setHighImpactPromiscuous() { mHighImpactPromiscuous = true; }
    public boolean isHighImpactPromiscuous() { return mHighImpactPromiscuous; }

    public boolean sameChromosome()
    {
        return mTranscripts[FS_UP].gene().chromosome().equals(mTranscripts[FS_DOWN].gene().chromosome());
    }

    public int distance()
    {
        if(!sameChromosome())
            return 0;

        return abs(mTranscripts[FS_UP].gene().position() - mTranscripts[FS_DOWN].gene().position());
    }

    public void setExonsSkipped(int exonsUp, int exonsDown)
    {
        mExonsSkipped[FS_UP] = exonsUp;
        mExonsSkipped[FS_DOWN] = exonsDown;
    }

    public int getExonsSkipped(boolean isUpstream) { return mExonsSkipped[fsIndex(isUpstream)]; }
    public int[] getExonsSkipped() { return mExonsSkipped; }

    public int getFusedExon(boolean isUpstream)
    {
        if(isExonic())
        {
            return isUpstream ? mTranscripts[FS_UP].ExonUpstream : mTranscripts[FS_DOWN].ExonDownstream;
        }

        if(isUpstream)
        {
            return max(mTranscripts[FS_UP].nextSpliceExonRank() - mExonsSkipped[FS_UP], 1);
        }
        else
        {
            return min(mTranscripts[FS_DOWN].nextSpliceExonRank() + mExonsSkipped[FS_DOWN], mTranscripts[FS_DOWN].exonCount());
        }
    }

    public int getBreakendExon(boolean isUpstream)
    {
        if(isExonic())
        {
            return isUpstream ? mTranscripts[FS_UP].ExonUpstream : mTranscripts[FS_DOWN].ExonDownstream;
        }

        return isUpstream ? mTranscripts[FS_UP].nextSpliceExonRank() : mTranscripts[FS_DOWN].nextSpliceExonRank();
    }

    public boolean isExonic()
    {
        return mTranscripts[FS_UP].isExonic() && mTranscripts[FS_DOWN].isExonic();
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

    public boolean nonDisruptiveChain()
    {
        if(mAnnotations != null && mAnnotations.chainInfo() != null)
            return mAnnotations.chainInfo().nonDisruptive();
        else
            return false;
    }

    public String svIdPair()
    {
        return String.format("%d_%d", mTranscripts[FS_UP].gene().id(), mTranscripts[FS_DOWN].gene().id());
    }

    public final String toString()
    {
        return String.format("%d: %s type(%s) reportable(%s reason=%s) phased(%s) SVs(%d & %d) up(%s:%d:%d exon=%d) down(%s:%d:%d exon=%d)",
                mId, name(), knownTypeStr(), mIsReportable, mReportableReason, mPhaseMatched,
                mTranscripts[FS_UP].gene().id(), mTranscripts[FS_DOWN].gene().id(),
                mTranscripts[FS_UP].gene().chromosome(), mTranscripts[FS_UP].gene().orientation(),
                mTranscripts[FS_UP].gene().position(), getFusedExon(true),
                mTranscripts[FS_DOWN].gene().chromosome(), mTranscripts[FS_DOWN].gene().orientation(),
                mTranscripts[FS_DOWN].gene().position(), getFusedExon(false));
    }
}
