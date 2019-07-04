package com.hartwig.hmftools.linx.cn;

import com.hartwig.hmftools.linx.types.SvBreakend;

public class HomLossEvent
{
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;
    public final String SegStart;
    public final String SegEnd;
    public final int StartSV; // the DB SvId
    public final int EndSV;

    private SvBreakend mBreakendStart;
    private SvBreakend mBreakendEnd;

    public HomLossEvent(
            final String chr,
            final long posStart,
            final long posEnd,
            final String segStart,
            final String segEnd,
            final int startSV,
            final int endSV)
    {
        Chromosome = chr;
        PosStart = posStart;
        PosEnd = posEnd;
        SegStart = segStart;
        SegEnd = segEnd;
        StartSV = startSV;
        EndSV = endSV;

        mBreakendStart = null;
        mBreakendEnd = null;
    }

    public void setBreakend(final SvBreakend breakend, boolean isStart)
    {
        if(isStart)
            mBreakendStart = breakend;
        else
            mBreakendEnd = breakend;
    }

    public void clearBreakends()
    {
        mBreakendStart = null;
        mBreakendEnd= null;
    }

    public final SvBreakend getBreakend(boolean isStart) { return isStart ? mBreakendStart : mBreakendEnd; }

    public boolean matchedBothSVs() { return mBreakendStart != null && mBreakendEnd != null; }
    public boolean sameSV() { return mBreakendStart != null && mBreakendStart.getSV() == mBreakendEnd.getSV(); }

    public boolean clustered()
    {
        return matchedBothSVs() && mBreakendStart.getCluster() == mBreakendEnd.getCluster();
    }

}
