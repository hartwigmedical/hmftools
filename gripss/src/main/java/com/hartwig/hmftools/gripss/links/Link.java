package com.hartwig.hmftools.gripss.links;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

public class Link
{
    public final String Id;
    private final Breakend[] mBreakends;

    public static final String LINK_TYPE_PAIR = "PAIR"; // used a single SV, not pairs of breakends from different SVs

    public Link(final String id, final Breakend first, final Breakend second)
    {
        Id = id;
        mBreakends = new Breakend[] { first, second };
    }

    public static Link from(final SvData sv)
    {
        return new Link(LINK_TYPE_PAIR, sv.breakendStart(), sv.breakendEnd());
    }

    public static Link from(final String linkId, final Breakend first, final Breakend second)
    {
        return new Link(linkId, first, second);
    }

    public String vcfId() { return mBreakends[SE_START].VcfId; }
    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public Breakend otherBreakend(final Breakend breakend) { return breakendStart() == breakend ?  breakendEnd() : breakendStart(); }

    public int minDistance()
    {
        if(breakendStart().sv().equals(breakendEnd().sv()))
            return breakendStart().sv().duplicationLength() + breakendStart().sv().insertSequenceLength();

        return abs(breakendStart().maxPosition() - breakendEnd().minPosition());
    }

    public int maxDistance()
    {
        if(breakendStart().sv().equals(breakendEnd().sv()))
            return breakendStart().sv().duplicationLength() + breakendStart().sv().insertSequenceLength();

        return abs(breakendStart().minPosition() - breakendEnd().maxPosition());
    }

    public Link reverse()
    {
        return new Link(Id, mBreakends[SE_END], mBreakends[SE_START]);
    }

    public String toString() { return String.format("%s<%s>%s distance(%d - %d)",
            breakendStart().VcfId, Id, breakendEnd().VcfId, minDistance(), maxDistance()); }
}
