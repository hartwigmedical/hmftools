package com.hartwig.hmftools.linx.cn;

import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;

public class HomLossEvent
{
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;
    public final String SegStart;
    public final String SegEnd;
    public final int StartSV; // the DB SvId
    public final int EndSV;

    private SvBreakend mBreakendStart;
    private SvBreakend mBreakendEnd;

    public HomLossEvent(
            final String chr,
            final int posStart,
            final int posEnd,
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
    public final SvCluster getCluster() { return matchedBothSVs() ? mBreakendStart.getCluster() : null; }

    public String toString()
    {
        return String.format("chr(%s) segs(%s -> %s) pos(%d -> %d) SVs(%s & %s)",
                Chromosome, SegStart, SegEnd, PosStart, PosEnd,
                mBreakendStart != null ? mBreakendStart.getSV().id() : "none", mBreakendEnd != null ? mBreakendEnd.getSV().id() : "none");
    }
}
