package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;

public class BaseDepth
{
    private final int[] mBaseRange;

    private int[] mDepth;
    private Map<Integer,Integer> mDepthMap;

    private static final int MIN_DEPTH_COUNT = 2;
    private static final int MAX_MAP_SIZE = 100000;
    private static final int DEFAULT_SIZE = 3000000; // largest observed gene collection range

    public BaseDepth()
    {
        mBaseRange = new int[SE_PAIR];
        mDepth = new int[DEFAULT_SIZE];
        mDepthMap = null;
    }

    public void initialise(final int[] baseRange)
    {
        mBaseRange[SE_START] = baseRange[SE_START];
        mBaseRange[SE_END] = baseRange[SE_END];

        if(length() > mDepth.length)
        {
            ISF_LOGGER.info("growing base depth: old({}) new({}) range({}-{})",
                    mDepth.length, length(), mBaseRange[SE_START], mBaseRange[SE_END]);

            mDepth = new int[length()];
        }
        else
        {
            // clear previous depth state
            for(int i = 0; i < mDepth.length; ++i)
                mDepth[i] = 0;
        }
    }

    public BaseDepth(final BaseDepth other, final Map<Integer,Integer> depthMap)
    {
        mBaseRange = new int[] { other.mBaseRange[SE_START], other.mBaseRange[SE_END] };
        mDepthMap = depthMap;
        mDepth = null;
    }

    public int length() { return mBaseRange[SE_END] - mBaseRange[SE_START] + 1; }

    public void processRead(final List<int[]> readCoords)
    {
        if(mDepth == null)
            return;

        for(final int[] readSection : readCoords)
        {
            int readStartPos = readSection[SE_START];
            int readEndPos = readSection[SE_END];

            if(readStartPos > mBaseRange[SE_END] || readEndPos < mBaseRange[SE_START])
                continue;

            // process this overlap
            int regionBaseIndex = readStartPos > mBaseRange[SE_START] ? readStartPos - mBaseRange[SE_START] : 0;
            int overlap = min(readEndPos, mBaseRange[SE_END]) - max(readStartPos, mBaseRange[SE_START]) + 1;

            if(regionBaseIndex + overlap > length())
            {
                ISF_LOGGER.error("baseDepth({}) read coords({} -> {}) regionBaseIndex({}) overlap({}) regionLength({})",
                        toString(), readStartPos, readEndPos, regionBaseIndex, overlap, length());
                return;
            }

            for(int j = regionBaseIndex; j < regionBaseIndex + overlap; ++j)
            {
                ++mDepth[j];
            }
        }
    }

    public Map<Integer,Integer> createPositionMap(final Set<Integer> candidateJunctions)
    {
        final Map<Integer,Integer> depthMap = Maps.newHashMap();

        for(final Integer position : candidateJunctions)
        {
            if(!positionWithin(position, mBaseRange[SE_START], mBaseRange[SE_END]))
                continue;

            int index = position - mBaseRange[SE_START];

            if(mDepth[index] >= MIN_DEPTH_COUNT)
                depthMap.put(position, mDepth[index]);
        }

        return depthMap;
    }

    public boolean hasPosition(int position) { return positionWithin(position, mBaseRange[SE_START], mBaseRange[SE_END]); }

    public int depthAtBase(int position)
    {
        if(!hasPosition(position))
            return 0;

        if(mDepthMap != null)
        {
            Integer depth = mDepthMap.get(position);
            return depth != null ? depth : 0;
        }

        int index = position - mBaseRange[SE_START];
        return mDepth[index];
    }

    public int basesWithDepth()
    {
        if(mDepth != null)
            return (int)Arrays.stream(mDepth).filter(x -> x >= MIN_DEPTH_COUNT).count();
        else
            return mDepthMap.size();
    }

    public double basesWithDepthPerc()
    {
        return basesWithDepth() / (double) mDepth.length;
    }

    public int maxDepth()
    {
        if(mDepth != null)
            return Arrays.stream(mDepth).max().orElse(0);
        else
            return mDepthMap.values().stream().mapToInt(x -> x).max().orElse(0);
    }

    public String toString()
    {
        return String.format("region(%d-%d)", mBaseRange[SE_START], mBaseRange[SE_END]);
    }

    public void collapse()
    {
        // move depth into a map if the coverage is relatively small over a large range
        mDepthMap = Maps.newHashMap();

        int index = 0;
        for(int pos = mBaseRange[SE_START]; pos <= mBaseRange[SE_END]; ++pos)
        {
            if(mDepth[index] >= MIN_DEPTH_COUNT)
                mDepthMap.put(pos, mDepth[index]);

            ++index;
        }

        if(mDepthMap.size() > MAX_MAP_SIZE)
        {
            ISF_LOGGER.warn("large map({} len={} perc={}) for baseDepth({})",
                    mDepthMap.size(), length(), String.format("%.2f", mDepthMap.size()/(double)length()), toString());
        }

        mDepth = null;
    }
}
