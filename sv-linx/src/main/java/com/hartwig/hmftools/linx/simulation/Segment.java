package com.hartwig.hmftools.linx.simulation;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.seIndex;

public class Segment
{
    public final int Id;

    private final boolean[] mHasLink;
    private final Segment[] mLinks;

    public Segment(int id, boolean hasStart, boolean hasEnd)
    {
        Id = id;
        mLinks = new Segment[SE_PAIR];
        mHasLink = new boolean[] {hasStart, hasEnd};
    }

    public boolean fullyLinked() { return !isLinkOpen(true) && !isLinkOpen(false); }

    public void clearLinks()
    {
        mLinks[SE_START] = null;
        mLinks[SE_END] = null;
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
