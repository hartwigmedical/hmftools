package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;

import java.util.Arrays;
import java.util.Map;

import com.google.common.collect.Maps;

public class BaseDepth
{
    public final int GeneCollection;
    public final String Chromosome;
    public final int[] BaseRange;

    private int[] mDepth;
    private Map<Integer,Integer> mDepthMap;

    private static final int MIN_DEPTH_COUNT = 2;
    private static final int RANGE_THRESHOLD = 3000; // little over the median gene length
    private static final int MAX_MAP_SIZE = 100000;

    public BaseDepth(final int geneCollection, final String chromosome, final int[] baseRange)
    {
        GeneCollection = geneCollection;
        Chromosome = chromosome;
        BaseRange = baseRange;
        mDepth = new int[BaseRange[SE_END] - BaseRange[SE_START] + 1];
        mDepthMap = null;
    }

    public void processRead(final ReadRecord read)
    {
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

    public void collapse()
    {
        // move depth into a map if the coverage is relatively small over a large range
        if(basesWithDepth() > MAX_MAP_SIZE)
            return;

        mDepthMap = Maps.newHashMap();
        int index = 0;
        for(int pos = BaseRange[SE_START]; pos <= BaseRange[SE_END]; ++pos)
        {
            if(mDepth[index] > MIN_DEPTH_COUNT)
                mDepthMap.put(pos, mDepth[index]);

            ++index;
        }

        mDepth = null;
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

    public int length() { return mDepth.length; }

    public int basesWithDepth()
    {
        return (int)Arrays.stream(mDepth).filter(x -> x > MIN_DEPTH_COUNT).count();
    }

    public double basesWithDepthPerc()
    {
        return basesWithDepth() / (double) mDepth.length;
    }

    public int maxDepth()
    {
        return Arrays.stream(mDepth).max().orElse(0);
    }

    public String toString()
    {
        return String.format("gc(%d) region(%s:%d-%d)", GeneCollection, Chromosome, BaseRange[SE_START], BaseRange[SE_END]);
    }


}
