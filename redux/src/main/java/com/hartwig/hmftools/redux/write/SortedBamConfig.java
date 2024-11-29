package com.hartwig.hmftools.redux.write;

public class SortedBamConfig
{
    // defaults
    private static final int READ_POS_CACHE_BUFFER = 30; // since reads with long soft-clips can lag the most recent duplicate position
    private static final int MIN_WRITE_COUNT = 10;
    private static final int POS_BUFFER_SIZE = 500;

    public final int PositionBuffer;
    public final int ReadPosCacheBuffer;
    public final int MinWriteCount;

    public SortedBamConfig()
    {
        this(POS_BUFFER_SIZE, READ_POS_CACHE_BUFFER, MIN_WRITE_COUNT);
    }

    public SortedBamConfig(final int positionBuffer, final int readPosCacheBuffer, final int minWriteCount)
    {
        PositionBuffer = positionBuffer;
        ReadPosCacheBuffer = readPosCacheBuffer;
        MinWriteCount = minWriteCount;
    }
}
