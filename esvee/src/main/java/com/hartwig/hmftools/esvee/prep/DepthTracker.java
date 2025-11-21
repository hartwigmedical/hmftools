package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.min;

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

    public double calcDepth(int position)
    {
        int index = windowIndex(position);
        return index != INVALID_INDEX && index < mWindowDepth.length ? mWindowDepth[index] / (double)mWindowSize : 0;
    }

    public void processRead(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        int indexStart = windowIndex(readStart);

        if(indexStart == INVALID_INDEX || indexStart >= mWindowDepth.length)
            return;

        int windowEnd = mRegion.start() + (indexStart + 1) * mWindowSize;
        int windowRange = min(readEnd, windowEnd) - readStart;

        mWindowDepth[indexStart] += windowRange;

        // assumes read spans at most 2 windows
        if(readEnd > windowEnd)
        {
            int indexEnd = indexStart + 1;

            if(indexEnd < mWindowDepth.length)
            {
                windowRange = readEnd - windowEnd;
                mWindowDepth[indexEnd] += windowRange;
            }
        }
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
