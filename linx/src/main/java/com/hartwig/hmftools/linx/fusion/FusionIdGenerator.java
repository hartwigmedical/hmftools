package com.hartwig.hmftools.linx.fusion;

public class FusionIdGenerator
{
    private int mNextFusionId;

    public FusionIdGenerator()
    {
        mNextFusionId = 0;
    }

    public void reset() { mNextFusionId = 0; }
    public int nextFusionId() { return mNextFusionId++; }
}
