package com.hartwig.hmftools.common.region;

import static java.lang.String.format;

public class UnmappingRegion extends BaseRegion
{
    private final int mMaxDepth;

    public UnmappingRegion(int posStart, int posEnd, int maxDepth)
    {
        super(posStart, posEnd);
        mMaxDepth = maxDepth;
    }

    public int maxDepth()
    {
        return mMaxDepth;
    }

    public String toString() { return format("%d-%d depth(%d)", start(), end(), mMaxDepth); }

    public static UnmappingRegion from(final BaseRegion region, int maxDepth)
    {
        return new UnmappingRegion(region.start(), region.end(), maxDepth);
    }

    public static UnmappingRegion from(final ChrBaseRegion region, int maxDepth)
    {
        return new UnmappingRegion(region.start(), region.end(), maxDepth);
    }
}
