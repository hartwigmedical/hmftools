package com.hartwig.hmftools.purple.somatic;

import htsjdk.variant.variantcontext.VariantContext;

public class KataegisWindow
{
    private final int mPosStart;
    private int mPosEnd;
    private int mCount;

    public KataegisWindow(final SomaticVariant variant)
    {
        mPosStart = variant.position();
        mPosEnd = mPosStart;
        mCount = 0;
    }

    public KataegisWindow(final KataegisWindow window)
    {
        mPosStart = window.start();
        mPosEnd = window.end();
        mCount = window.mCount;
    }

    public void add(final SomaticVariant variant)
    {
        mCount++;
        mPosEnd = variant.position();
    }

    public int count()
    {
        return mCount;
    }

    public int start()
    {
        return mPosEnd;
    }
    public int end()
    {
        return mPosEnd;
    }

    boolean isViable(int minCount, long maxAverageDistance)
    {
        return count() >= minCount && averageDistance() <= maxAverageDistance;
    }

    long averageDistance()
    {
        return mCount == 0 || mCount == 1 ? 0 : Math.round(1d * (mPosEnd - mPosStart) / (mCount - 1));
    }
}
