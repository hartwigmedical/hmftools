package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
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
    private final int[] mWindowReadCounts;

    public DepthTracker(final BaseRegion region, final int windowSize)
    {
        mRegion = region;
        mWindowSize = windowSize;

        int regionLength = mRegion.baseLength();
        int windows = (int)ceil(regionLength / mWindowSize);

        mWindowDepth = new long[windows];
        mWindowReadCounts = new int[windows];
    }

    public double calcDepth(int position)
    {
        int index = windowIndex(position);
        return index != INVALID_INDEX && index < mWindowDepth.length ? mWindowDepth[index] / (double)mWindowSize : 0;
    }

    public int windowReadCount(int position)
    {
        int index = windowIndex(position);
        return index != INVALID_INDEX && index < mWindowReadCounts.length ? mWindowReadCounts[index] : 0;
    }

    public void processRead(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
            return;

        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        int indexStart = windowIndex(readStart);

        if(indexStart == INVALID_INDEX || indexStart >= mWindowDepth.length)
            return;

        int windowEnd = mRegion.start() + (indexStart + 1) * mWindowSize;
        int windowRange = min(readEnd, windowEnd) - readStart;

        if(windowRange < 0 || windowRange > mWindowSize)
        {
            SV_LOGGER.error("depth range failed: read({}:{}-{}) window(idx={} end={}) range({})",
                    readStart, readEnd, indexStart, windowEnd, windowRange);
            return;
        }

        mWindowDepth[indexStart] += windowRange;
        ++mWindowReadCounts[indexStart];

        // assumes read spans at most 2 windows
        if(readEnd > windowEnd)
        {
            int indexEnd = indexStart + 1;

            if(indexEnd < mWindowDepth.length)
            {
                windowRange = readEnd - windowEnd;

                if(windowRange < 0 || windowRange > mWindowSize)
                {
                    SV_LOGGER.error("depth range failed: read({}:{}-{}) window(idx={} end={}) range({})",
                            readStart, readEnd, indexEnd, windowEnd, windowRange);
                    return;
                }

                mWindowDepth[indexEnd] += windowRange;
                ++mWindowReadCounts[indexEnd];
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
