package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import java.util.Map;

import com.google.common.collect.Maps;

public class RegionAllocator
{
    private final int mBlockSize;
    private final int mBlockArea;

    private Map<String,Boolean> mAllocations;

    public static int DEFAULT_REGION_GRID_SIZE = 100;

    public RegionAllocator(int blockSize)
    {
        mBlockSize = blockSize;
        mBlockArea= blockSize * blockSize;

        mAllocations = Maps.newHashMap();
    }

    public int blockSize() { return mBlockSize; }
    public int allocationCount() { return mAllocations.size(); }

    public int baseToIndex(long base)
    {
        double rawIndex = base/(double)mBlockSize;
        return (int)floor(rawIndex);
    }

    private static final String DELIM = "_";

    public static String basePairToKey(int region1, int region2)
    {
        return String.valueOf(region1) + DELIM + String.valueOf(region2);
    }

    public boolean canAllocate(int region1, int region2)
    {
        final String key = basePairToKey(region1, region2);
        Boolean allocated = mAllocations.get(key);
        if(allocated != null)
            return false;

        mAllocations.put(key, true);
        return true;
    }

    public long allocateBases(
            long lowerStart, long lowerEnd, long upperStart, long upperEnd,
            long minBucketLen, long maxBucketLen, boolean requireAllocation)
    {
        if(lowerStart >= lowerEnd || upperStart >= upperEnd)
            return 0;

        int overlapCount = 0;

        for (long base = lowerStart; base <= lowerEnd;)
        {
            int baseIndex = baseToIndex(base);

            int upperStartIndex = baseToIndex(max(upperStart, base + minBucketLen));
            int upperEndIndex = baseToIndex(min(upperEnd, base + maxBucketLen));

            if(upperStartIndex <= upperEndIndex)
            {
                for (int upperIndex = upperStartIndex; upperIndex <= upperEndIndex; ++upperIndex)
                {
                    if (!requireAllocation || canAllocate(baseIndex, upperIndex))
                    {
                        ++overlapCount;
                    }
                }
            }

            base += mBlockSize;
        }

        return overlapCount * mBlockArea;
    }

    public long allocateBases(long lowerStart, long lowerEnd, long upperStart, long upperEnd)
    {
        if(lowerStart >= lowerEnd || upperStart >= upperEnd)
            return 0;

        int overlapCount = 0;

        for (long base = lowerStart; base <= lowerEnd;)
        {
            int baseIndex = baseToIndex(base);

            int upperStartIndex = baseToIndex(upperStart);
            int upperEndIndex = baseToIndex(upperEnd);

            if(upperStartIndex <= upperEndIndex)
            {
                for (int upperIndex = upperStartIndex; upperIndex <= upperEndIndex; ++upperIndex)
                {
                    if (canAllocate(baseIndex, upperIndex))
                    {
                        ++overlapCount;
                    }
                }
            }

            base += mBlockSize;
        }

        return overlapCount * mBlockArea;
    }

    public void reset()
    {
        mAllocations.clear();
    }


}
