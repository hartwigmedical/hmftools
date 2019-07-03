package com.hartwig.hmftools.linx.cn;

import com.hartwig.hmftools.linx.types.SvBreakend;

public class HomLossEvent
{
    public final String SampleId;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;
    public final String SegStart;
    public final String SegEnd;
    public final int StartSV; // the SV data's ID
    public final int EndSV;

    private SvBreakend mBreakendStart;
    private SvBreakend mBreakendEnd;

    public HomLossEvent(
            final String sampleId,
            final String chr,
            final long posStart,
            final long posEnd,
            final String segStart,
            final String segEnd,
            final int startSV,
            final int endSV)
    {
        SampleId = sampleId;
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

    public final SvBreakend getBreakend(boolean isStart) { return isStart ? mBreakendStart : mBreakendEnd; }

    public boolean matchedBothSVs() { return mBreakendStart != null && mBreakendEnd != null; }

}
