package com.hartwig.hmftools.bamtools.remapper.testutilities;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class BaseRegion extends ChrBaseRegion
{
    public final byte[] mBases;
    public final int mChromosomeIndex;

    public BaseRegion(final int chromosomeIndex, int start, int end, byte[] bases)
    {
        super("chr" + (chromosomeIndex + 1), start, end);
        mBases = bases;
        mChromosomeIndex = chromosomeIndex;
    }
}
