package com.hartwig.hmftools.cobalt.lowcov;

public class LowCovBucket
{
    public final int startPosition;
    public final int endPosition;

    public final int bucketPosition;

    public LowCovBucket(int startPosition, int endPosition, int midPosition)
    {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.bucketPosition = midPosition;
    }
}
