package com.hartwig.hmftools.svtools.simulation;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.seIndex;

public class Segment
{
    public final int Id;
    public final boolean EndSegment; // either the start or end of the chain, always retained and only linked on one end

    private final boolean[] mHasLink;
    private final Segment[] mLinks;

    public Segment(int id, boolean hasStart, boolean hasEnd)
    {
        Id = id;
        mLinks = new Segment[SE_PAIR];
        mHasLink = new boolean[] {hasStart, hasEnd};

        EndSegment = !hasStart || !hasEnd;
    }

    public static Segment from(final Segment other)
    {
        return new Segment(other.Id, other.hasLink(true), other.hasLink(false));
    }

    public boolean fullyLinked() { return !isLinkOpen(true) && !isLinkOpen(false); }

    public void clearLinks()
    {
        mLinks[SE_START] = null;
        mLinks[SE_END] = null;
    }

    public boolean hasLink(boolean isStart)
    {
        return mHasLink[seIndex(isStart)];
    }
    public boolean isLinkOpen(boolean isStart)
    {
        return mHasLink[seIndex(isStart)] && mLinks[seIndex(isStart)] == null;
    }

    public Segment getLink(boolean isStart) { return mLinks[seIndex(isStart)]; }

    public void setLink(final Segment link, boolean isStart)
    {
        mLinks[seIndex(isStart)] = link;
    }

}
