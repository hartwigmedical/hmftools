package com.hartwig.hmftools.esvee.vcfcompare.line;

import com.hartwig.hmftools.esvee.vcfcompare.VariantBreakend;

public class LineLink
{
    public final VariantBreakend mPolyASite;
    public final VariantBreakend mOtherSite;
    public final LineLinkType mType;
    
    public LineLink(VariantBreakend polyASite, VariantBreakend otherSite, LineLinkType type)
    {
        checkPolyASite(polyASite);

        mPolyASite = polyASite;
        mOtherSite = otherSite;
        mType = type;
    }

    private static void checkPolyASite(VariantBreakend polyASite)
    {
        if(!polyASite.hasPolyATail())
        {
            throw new IllegalStateException(String.format("breakend(%s) does not have a poly A tail. Insert sequence: %s",
                            polyASite, polyASite.InsertSequence));
        }
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
