package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMRecord;

public class DepthTracker
{
    private final BaseRegion mRegion;
    private final int mWindowSize;
    private final long[] mWindowDepth;

    public DepthTracker(final BaseRegion region, final int windowSize)
    {
        mRegion = region;
        mWindowSize = windowSize;

        int regionLength = mRegion.baseLength();
        int windows = (int)ceil(regionLength / mWindowSize);

        mWindowDepth = new long[windows];
    }

    public double calcDepth(int position)
    {
        int index = windowIndex(position);
        return index != INVALID_INDEX && index < mWindowDepth.length ? mWindowDepth[index] / (double)mWindowSize : 0;
    }

    public void processRead(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
            return;

        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        int windowIndex = windowIndex(readStart);

        if(windowIndex == INVALID_INDEX || windowIndex >= mWindowDepth.length)
            return;

        // apply the read's bases across all applicable windows
        // int windowStart = readStart;
        int windowStart = mRegion.start() + windowIndex * mWindowSize;;
        int windowEnd = windowStart + mWindowSize - 1;

        while(true)
        {
            int windowRange = min(windowEnd, readEnd) - max(readStart, windowStart) + 1;

            if(windowRange < 0 || windowRange > mWindowSize)
            {
                SV_LOGGER.error("depth range failed: read({}-{}) window(idx={} {}-{}) read range({})",
                        readStart, readEnd, windowIndex, windowStart, windowEnd, windowRange);
                return;
            }

            mWindowDepth[windowIndex] += windowRange;

            if(readEnd <= windowEnd)
                break;

            ++windowIndex;

            if(windowIndex >= mWindowDepth.length)
                break;

            windowStart = windowEnd + 1;
            windowEnd = windowStart + mWindowSize - 1;
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
