package com.hartwig.hmftools.linx.simulation;

public class Segment
{
    public final int Id;
    public final boolean HasStart;
    public final boolean HasEnd;

    private Segment mStartLink;
    private Segment mEndLink;

    public Segment(int id, boolean hasStart, boolean hasEnd)
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

    public Segment getLink(boolean isStart) { return isStart ? mStartLink : mEndLink; }

    public void setLink(final Segment link, boolean isStart)
    {
        if(isStart)
            mStartLink = link;
        else
            mEndLink = link;
    }

}
