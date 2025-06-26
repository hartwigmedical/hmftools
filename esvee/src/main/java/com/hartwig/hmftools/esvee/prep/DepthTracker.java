package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class DepthTracker
{
    private final BaseRegion mRegion;
    private final int mWindowSize;
    private final int[] mWindowDepth;

    public DepthTracker(final BaseRegion region, final int windowSize)
    {
        mRegion = region;
        mWindowSize = windowSize;

        int regionLength = mRegion.baseLength();
        int windows = (int)ceil(regionLength / mWindowSize);

        mWindowDepth = new int[windows];
    }

    public int calcDepth(int position)
    {
        int index = windowIndex(position);

        return index != INVALID_INDEX && index < mWindowDepth.length ? mWindowDepth[index] : 0;
    }

    public void processRead(final SAMRecord read)
    {
        int indexStart = windowIndex(read.getAlignmentStart());

        if(indexStart == INVALID_INDEX || indexStart >= mWindowDepth.length)
            return;

        ++mWindowDepth[indexStart];

        int indexEnd = windowIndex(read.getAlignmentEnd());

        if(indexEnd > indexStart && indexEnd < mWindowDepth.length)
            ++mWindowDepth[indexEnd];
    }

    private static final int INVALID_INDEX = -1;

    private int windowIndex(int position)
    {
        if(!positionWithin(position, mRegion.start(), mRegion.end()))
            return INVALID_INDEX;

        int positionOffset = position - mRegion.start();

        return (int)floor(positionOffset / mWindowSize);
    }
}
