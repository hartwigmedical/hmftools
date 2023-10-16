package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_POS_BUFFER_SIZE;

public class SortedBamConfig
{
    // defaults
    private static final double CAPACITY_CHECK_PERCENT = 0.01;
    private static final double CAPACITY_GROW_PERCENT = 0.9;
    private static final double CAPACITY_SHRINK_PERCENT = 0.5;

    private static final int READ_POS_CACHE_BUFFER = 30; // since reads with long soft-clips can lag the most recent duplicate position
    private static final int MIN_SORT_WRITE_COUNT = 10;
    private static final int BAM_READ_CACHE_BUFFER = 20000;

    public final int Capacity;
    public final int PositionBuffer;
    public final double CapacityCheckPercent;
    public final double CapacityGrowPercent;
    public final double CapacityShrinkPercent;
    public final int ReadPosCacheBuffer;
    public final int MinSortWriteCount;

    public SortedBamConfig()
    {
        this(BAM_READ_CACHE_BUFFER, DEFAULT_POS_BUFFER_SIZE, READ_POS_CACHE_BUFFER, MIN_SORT_WRITE_COUNT,
                CAPACITY_CHECK_PERCENT, CAPACITY_GROW_PERCENT, CAPACITY_SHRINK_PERCENT);
    }

    public SortedBamConfig(
            final int capacity, final int positionBuffer, final int readPosCacheBuffer, final int minSortWriteCount,
            final double capacityCheckPercent, final double capacityGrowPercent, final double capacityShrinkPercent)
    {
        Capacity = capacity;
        PositionBuffer = positionBuffer;
        CapacityCheckPercent = capacityCheckPercent;
        CapacityGrowPercent = capacityGrowPercent;
        CapacityShrinkPercent = capacityShrinkPercent;
        ReadPosCacheBuffer = readPosCacheBuffer;
        MinSortWriteCount = minSortWriteCount;
    }
}
