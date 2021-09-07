package com.hartwig.hmftools.linx.cn;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.Q_ARM;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class LohEvent
{
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;
    public final String SegStart;
    public final String SegEnd;
    public final int SegCount;
    public final int StartSV; // the DB SvId
    public final int EndSV;

    private final List<HomLossEvent> mHomLossEvents;

    private SvCNData mStartCnData;
    private SvCNData mEndCnData;

    private boolean mIsValid;
    private SvBreakend mBreakendStart;
    private SvBreakend mBreakendEnd;

    public static final int CN_DATA_NO_SV = -1;

    public LohEvent(
            final String chr,
            final int posStart,
            final int posEnd,
            final String segStart,
            final String segEnd,
            final int segCount,
            final int startSV,
            final int endSV)
    {
        Chromosome = chr;
        PosStart = posStart;
        PosEnd = posEnd;
        SegStart = segStart;
        SegEnd = segEnd;
        SegCount = segCount;
        StartSV = startSV;
        EndSV = endSV;

        mStartCnData = null;
        mEndCnData = null;
        mBreakendStart = null;
        mBreakendEnd = null;
        mHomLossEvents = Lists.newArrayList();
        mIsValid = true;
    }

    public void setBreakend(final SvBreakend breakend, boolean isStart)
    {
        if(isStart)
            mBreakendStart = breakend;
        else
            mBreakendEnd = breakend;
    }

    public final SvBreakend getBreakend(boolean isStart) { return isStart ? mBreakendStart : mBreakendEnd; }

    public void setCnData(final SvCNData startCnData, final SvCNData endCnData)
    {
        mStartCnData = startCnData;
        mEndCnData = endCnData;
    }

    public final SvCNData getCnData(boolean isStart) { return isStart ? mStartCnData : mEndCnData; }

    public boolean matchedBothSVs() { return mBreakendStart != null && mBreakendEnd != null; }
    public boolean sameSV() { return mBreakendStart != null && mBreakendStart.getSV() == mBreakendEnd.getSV(); }
    public int length() { return PosEnd - PosStart; }

    public void addHomLossEvents(final List<HomLossEvent> events)
    {
        if(events == null)
            return;

        mHomLossEvents.addAll(events);
        mIsValid = !hasIncompleteHomLossEvents();
    }

    public final List<HomLossEvent> getHomLossEvents() { return mHomLossEvents; }

    public boolean hasIncompleteHomLossEvents()
    {
        return mHomLossEvents.stream().anyMatch(x -> !x.matchedBothSVs() || !x.sameSV());
    }

    public boolean clustered()
    {
        return matchedBothSVs() && mBreakendStart.getCluster() == mBreakendEnd.getCluster();
    }

    public boolean isValid() { return mIsValid; }
    public void setIsValid(boolean toggle) { mIsValid = toggle; }

    public boolean armLoss()
    {
        return armLoss(P_ARM) || armLoss(Q_ARM);
    }

    public boolean armLoss(final ChromosomeArm arm)
    {
        if(arm == P_ARM && SegStart.equals(TELOMERE.toString()) && SegEnd.equals(CENTROMERE.toString()))
            return true;

        if(arm == Q_ARM && SegStart.equals(CENTROMERE.toString()) && SegEnd.equals(TELOMERE.toString()))
            return true;

        return false;
    }

    public boolean telomereLoss()
    {
        return SegStart.equals(TELOMERE.toString()) || SegEnd.equals(TELOMERE.toString());
    }

    public boolean centromereLoss()
    {
        return SegStart.equals(CENTROMERE.toString()) || SegEnd.equals(CENTROMERE.toString());
    }

    public boolean chromosomeLoss()
    {
        return SegStart.equals(TELOMERE.toString()) && SegEnd.equals(TELOMERE.toString());
    }

    public boolean wholeArmLoss() { return armLoss() || chromosomeLoss(); }

    public boolean isSvEvent(boolean useStart)
    {
        return useStart ? StartSV != CN_DATA_NO_SV : EndSV != CN_DATA_NO_SV;
    }
    public boolean isSvEvent() { return isSvEvent(true) || isSvEvent(false); }
    public boolean doubleSvEvent() { return isSvEvent(true) && isSvEvent(false); }

    public String toString()
    {
        return String.format("chr(%s) segs(%s -> %s) pos(%d -> %d) SVs(%s & %s)",
                Chromosome, SegStart, SegEnd, PosStart, PosEnd,
                mBreakendStart != null ? mBreakendStart.getSV().id() : "none", mBreakendEnd != null ? mBreakendEnd.getSV().id() : "none");
    }

}
