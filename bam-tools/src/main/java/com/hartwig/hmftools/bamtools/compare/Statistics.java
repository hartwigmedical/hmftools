package com.hartwig.hmftools.bamtools.compare;

public class Statistics
{
    public int RefReadCount;
    public int NewReadCount;
    public int DiffCount;

    public Statistics()
    {
        RefReadCount = 0;
        NewReadCount = 0;
        DiffCount = 0;
    }

    public void merge(final Statistics other)
    {
        RefReadCount += other.RefReadCount;
        NewReadCount += other.NewReadCount;
        DiffCount += other.DiffCount;
    }
}
