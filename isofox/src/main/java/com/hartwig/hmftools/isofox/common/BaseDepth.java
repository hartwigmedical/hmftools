package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;

public class BaseDepth
{
    public final int GeneCollection;
    public final String Chromosome;
    public final int[] BaseRange;

    private int[] mDepth;
    private Map<Integer,Integer> mDepthMap;

    private static final int MIN_DEPTH_COUNT = 2;
    private static final int MAX_MAP_SIZE = 100000;

    public BaseDepth(final int geneCollection, final String chromosome, final int[] baseRange)
    {
        GeneCollection = geneCollection;
        Chromosome = chromosome;
        BaseRange = baseRange;
        mDepth = new int[BaseRange[SE_END] - BaseRange[SE_START] + 1];
        mDepthMap = null;
    }

    public BaseDepth(final BaseDepth other, final Map<Integer,Integer> depthMap)
    {
        GeneCollection = other.GeneCollection;
        Chromosome = other.Chromosome;
        BaseRange = other.BaseRange;
        mDepthMap = depthMap;
        mDepth = null;
    }

    public void processRead(final ReadRecord read)
    {
        if(mDepth == null)
            return;

        for(final int[] readSection : read.getMappedRegionCoords())
        {
            int readStartPos = readSection[SE_START];
            int readEndPos = readSection[SE_END];

            if (readStartPos > BaseRange[SE_END] || readEndPos < BaseRange[SE_START])
                continue;

            // process this overlap
            int regionBaseIndex = readStartPos > BaseRange[SE_START] ? readStartPos - BaseRange[SE_START] : 0;
            int overlap = min(readEndPos, BaseRange[SE_END]) - max(readStartPos, BaseRange[SE_START]) + 1;

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

    public void clearDepth()
    {
        mDepth = null;
    }

    public Map<Integer,Integer> createPositionMap(final Set<Integer> candidateJunctions)
    {
        final Map<Integer,Integer> depthMap = Maps.newHashMap();

        for(final Integer position : candidateJunctions)
        {
            if(!positionWithin(position, BaseRange[SE_START], BaseRange[SE_END]))
                continue;

            int index = position - BaseRange[SE_START];

            if(mDepth[index] > MIN_DEPTH_COUNT)
                depthMap.put(position, mDepth[index]);
        }

        return depthMap;
    }

    public int depthAtBase(int position)
    {
        if(!positionWithin(position, BaseRange[SE_START], BaseRange[SE_END]))
            return 0;

        if(mDepthMap != null)
        {
            Integer depth = mDepthMap.get(position);
            return depth != null ? depth : 0;
        }

        int index = position - BaseRange[SE_START];
        return mDepth[index];
    }

    public int length() { return BaseRange[SE_END] - BaseRange[SE_START] + 1; }

    public int basesWithDepth()
    {
        if(mDepth != null)
            return (int)Arrays.stream(mDepth).filter(x -> x > MIN_DEPTH_COUNT).count();
        else
            return (int)mDepthMap.values().stream().filter(x -> x > MIN_DEPTH_COUNT).count();
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
        return String.format("gc(%d) region(%s:%d-%d)", GeneCollection, Chromosome, BaseRange[SE_START], BaseRange[SE_END]);
    }

    public void collapse()
    {
        // move depth into a map if the coverage is relatively small over a large range
        mDepthMap = Maps.newHashMap();

        int index = 0;
        for(int pos = BaseRange[SE_START]; pos <= BaseRange[SE_END]; ++pos)
        {
            if(mDepth[index] > MIN_DEPTH_COUNT)
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
