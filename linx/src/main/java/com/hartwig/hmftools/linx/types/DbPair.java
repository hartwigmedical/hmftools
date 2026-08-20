package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.common.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;

public class DbPair
{
    private final SvBreakend mLowerBreakend;
    private final SvBreakend mUpperBreakend;
    private final int mLinkLength;

    public DbPair(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        if(breakend1.position() < breakend2.position()
        || (breakend1.position() == breakend2.position() && breakend1.orientation() == ORIENT_FWD))
        {
            mLowerBreakend = breakend1;
            mUpperBreakend = breakend2;
        }
        else
        {
            mLowerBreakend = breakend2;
            mUpperBreakend = breakend1;
        }

        // exact base is considered a 1 base overlap, 1 base apart a zero-length DB
        int length = mUpperBreakend.position() - mLowerBreakend.position() - 1;

        if(mLowerBreakend.orientation() == ORIENT_FWD)
        {
            mLinkLength = length;
        }
        else
        {
            mLinkLength = -length;
        }
    }

    public SvBreakend lower() { return mLowerBreakend; }
    public SvBreakend upper() { return mUpperBreakend; }

    public SvVarData lowerSV() { return mLowerBreakend.getSV(); }
    public SvVarData upperSV() { return mUpperBreakend.getSV(); }

    public boolean lowerLinkOnStart() { return mLowerBreakend.usesStart(); }
    public boolean upperLinkOnStart() { return mUpperBreakend.usesStart(); }

    public SvBreakend getBreakend(int se) { return getBreakend(isStart(se)); }

    public SvBreakend getBreakend(boolean isStart) { return isStart ? mLowerBreakend : mUpperBreakend; }

    public String chromosome() { return mLowerBreakend.chromosome(); }

    public int length() { return mLinkLength; }

    public boolean hasBreakend(final SvVarData var, boolean useStart)
    {
        return (var == mLowerBreakend.getSV() && mLowerBreakend.usesStart() == useStart)
            || (var == mUpperBreakend.getSV() && mUpperBreakend.usesStart() == useStart);
    }

    public boolean hasBreakend(final SvBreakend breakend)
    {
        return breakend == mLowerBreakend || breakend == mUpperBreakend;
    }

    public String toString()
    {
        return String.format("%s %s:%d:%s & %s %s:%d:%s",
                mLowerBreakend.getSV().id(), mLowerBreakend.chromosome(), mLowerBreakend.position(),
                mLowerBreakend.usesStart() ? "start" : "end",
                mUpperBreakend.getSV().id(), mUpperBreakend.chromosome(), mUpperBreakend.position(),
                mUpperBreakend.usesStart() ? "start" : "end");
    }

    public SvVarData getOtherSV(final SvVarData var) { return mLowerBreakend.getSV() == var ? mUpperBreakend.getSV() : mLowerBreakend.getSV(); }

    public SvBreakend getOtherBreakend(final SvBreakend breakend)
    {
        return mLowerBreakend == breakend ? mUpperBreakend : mLowerBreakend;
    }
}
