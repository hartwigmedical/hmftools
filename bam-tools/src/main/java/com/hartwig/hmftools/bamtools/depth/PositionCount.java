package com.hartwig.hmftools.bamtools.depth;

class PositionCount
{
    public int Position;
    public int Count;
    public int DepthMin;
    public int DepthMax;

    public PositionCount(final int position, int depthMin, int depthMax)
    {
        Position = position;
        Count = 1;
        DepthMax = depthMax;
        DepthMin = depthMin;
    }
}
