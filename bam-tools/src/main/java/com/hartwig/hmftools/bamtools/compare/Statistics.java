package com.hartwig.hmftools.bamtools.compare;

public class Statistics
{
    public long OrigReadCount;
    public long NewReadCount;
    public long Matched;
    public long DiffCount;
    public long ValueDiffCount;
    public long OrigOnlyCount;
    public long NewOnlyCount;

    public void recordOrigOnly()
    {
        ++OrigOnlyCount;
        ++DiffCount;
    }

    public void recordNewOnly()
    {
        ++NewOnlyCount;
        ++DiffCount;
    }

    public void recordValue()
    {
        ++ValueDiffCount;
        ++DiffCount;
    }

    public void merge(final Statistics other)
    {
        OrigReadCount += other.OrigReadCount;
        NewReadCount += other.NewReadCount;
        Matched += other.Matched;
        DiffCount += other.DiffCount;
        ValueDiffCount += other.ValueDiffCount;
        OrigOnlyCount += other.OrigOnlyCount;
        NewOnlyCount += other.NewOnlyCount;
    }
}
