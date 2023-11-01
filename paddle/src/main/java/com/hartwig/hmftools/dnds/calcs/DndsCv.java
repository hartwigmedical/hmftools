package com.hartwig.hmftools.dnds.calcs;

import static java.lang.Math.max;

public class DndsCv
{
    private final int mCount; // was n
    public final double mWCv; // was wCv

    public DndsCv(final int count, final double wCv)
    {
        mCount = count;
        mWCv = wCv;
    }

    public double expectedDrivers(int variantCount)
    {
        double probability = mCount > 0 & mWCv > 0 ? max((mWCv - 1) / mWCv, 0) : 0;
        return variantCount * probability;

    }
}