package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class HighDepthRegion extends BaseRegion
{
    private final int mMaxDepth;

    public HighDepthRegion(int posStart, int posEnd, int maxDepth)
    {
        super(posStart, posEnd);
        mMaxDepth = maxDepth;
    }

    public int maxDepth()
    {
        return mMaxDepth;
    }

    public String toString() { return format("%d-%d depth(%d)", start(), end(), mMaxDepth); }

    public static HighDepthRegion from(final BaseRegion region, int maxDepth)
    {
        return new HighDepthRegion(region.start(), region.end(), maxDepth);
    }

    public static HighDepthRegion from(final ChrBaseRegion region, int maxDepth)
    {
        return new HighDepthRegion(region.start(), region.end(), maxDepth);
    }
}
