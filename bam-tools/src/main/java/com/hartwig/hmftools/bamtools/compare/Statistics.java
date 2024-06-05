package com.hartwig.hmftools.bamtools.compare;

public class Statistics
{
    public long OrigReadCount;
    public long NewReadCount;
    public long DiffCount;

    public Statistics()
    {
        OrigReadCount = 0;
        NewReadCount = 0;
        DiffCount = 0;
    }

    public void merge(final Statistics other)
    {
        OrigReadCount += other.OrigReadCount;
        NewReadCount += other.NewReadCount;
        DiffCount += other.DiffCount;
    }
}
