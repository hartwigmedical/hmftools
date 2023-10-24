package com.hartwig.hmftools.markdups.common;

import com.hartwig.hmftools.common.region.BaseRegion;

public class HighDepthRegion extends BaseRegion
{
    private final int mMaxDepth;

    public HighDepthRegion(final int posStart, final int posEnd, final int maxDepth)
    {
        super(posStart, posEnd);
        mMaxDepth = maxDepth;
    }

    public int maxDepth()
    {
        return mMaxDepth;
    }
}
