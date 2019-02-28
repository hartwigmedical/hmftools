package com.hartwig.hmftools.svanalysis.simulation;

public class SvSimSegment
{
    public final int Id;
    public final boolean HasStart;
    public final boolean HasEnd;

    private SvSimSegment mStartLink;
    private SvSimSegment mEndLink;

    public SvSimSegment(int id, boolean hasStart, boolean hasEnd)
    {
        Id = id;
        mStartLink = null;
        mEndLink = null;
        HasStart = hasStart;
        HasEnd = hasEnd;
    }

    public boolean fullyLinked() { return !isLinkOpen(true) && !isLinkOpen(false); }
    public void clearLinks()
    {
        mStartLink = null;
        mEndLink = null;
    }

    public boolean isLinkOpen(boolean isStart)
    {
        return isStart ? HasStart && mStartLink == null: HasEnd && mEndLink == null;
    }

    public SvSimSegment getLink(boolean isStart) { return isStart ? mStartLink : mEndLink; }

    public void setLink(final SvSimSegment link, boolean isStart)
    {
        if(isStart)
            mStartLink = link;
        else
            mEndLink = link;
    }

}
