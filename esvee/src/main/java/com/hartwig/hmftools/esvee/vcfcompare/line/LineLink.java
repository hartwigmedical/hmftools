package com.hartwig.hmftools.esvee.vcfcompare.line;

import static java.lang.String.format;

import com.hartwig.hmftools.esvee.vcfcompare.Breakend;

public class LineLink
{
    public final Breakend mPolyASite;
    public final Breakend mOtherSite;
    public final LineLinkType mType;
    
    public LineLink(Breakend polyASite, Breakend otherSite, LineLinkType type)
    {
        checkPolyASite(polyASite);

        mPolyASite = polyASite;
        mOtherSite = otherSite;
        mType = type;
    }

    private static void checkPolyASite(Breakend polyASite)
    {
        if(!polyASite.hasPolyATail())
        {
            throw new IllegalStateException(format("breakend(%s) does not have a poly A tail. Insert sequence: %s",
                            polyASite, polyASite.InsertSequence));
        }
    }

    public boolean hasRemotes()
    {
        return !mPolyASite.sv().id().equals(mOtherSite.sv().id());
    }

    public boolean polyAHasRemote()
    {
        return !mPolyASite.isSgl() && hasRemotes();
    }

    public boolean otherHasRemote()
    {
        return !mOtherSite.isSgl() && hasRemotes();
    }

    @Override
    public String toString()
    {
        return format("polyASite(%s %s) otherSite(%s %s)",
                mPolyASite.sv().id(), mPolyASite.coordStr(), mOtherSite.sv().id(), mOtherSite.coordStr());
    }
}
