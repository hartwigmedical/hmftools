package com.hartwig.hmftools.cobalt.lowcov;

public class LowCovBucket
{
    public final int StartPosition;
    public final int EndPosition;
    public final int BucketPosition;

    public LowCovBucket(int startPosition, int endPosition, int midPosition)
    {
        StartPosition = startPosition;
        EndPosition = endPosition;
        BucketPosition = midPosition;
    }

    @Override
    public String toString()
    {
        return "LowCovBucket{" +
                "StartPosition=" + StartPosition +
                ", EndPosition=" + EndPosition +
                ", BucketPosition=" + BucketPosition +
                '}';
    }
}
