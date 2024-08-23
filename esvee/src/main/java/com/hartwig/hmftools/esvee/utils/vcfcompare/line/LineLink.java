package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;

public class LineLink
{
    VariantBreakend mPolyASite;
    VariantBreakend mOtherSite;

    public LineLink(VariantBreakend polyASite, VariantBreakend otherSite)
    {
        mPolyASite = polyASite;
        mOtherSite = otherSite;
    }

    public boolean hasRemotes()
    {
        return !mPolyASite.eventId().equals(mOtherSite.eventId());
    }

    public boolean polyAHasRemote()
    {
        return !mPolyASite.isSingle() && hasRemotes();
    }

    public boolean otherHasRemote()
    {
        return !mOtherSite.isSingle() && hasRemotes();
    }

    @Override
    public String toString()
    {
        return String.format("polyASite(%s %s) otherSite(%s %s)", mPolyASite.Id, mPolyASite.coordStr(), mOtherSite.Id, mOtherSite.coordStr());
    }
}
