package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

public class EvidenceStats
{
    public long ReadCount;
    public long NoVariantReadCount;
    public final long[] SupportCounts;

    public int PartitionCount;
    public int SliceCount;
    public int SliceLength;

    public EvidenceStats()
    {
        ReadCount = 0;
        NoVariantReadCount = 0;
        PartitionCount = 0;
        SliceCount = 0;
        SliceLength = 0;
        SupportCounts = new long[ReadMatchType.values().length];
    }

    public void merge(final EvidenceStats other)
    {
        ReadCount += other.ReadCount;
        NoVariantReadCount += other.NoVariantReadCount;
        PartitionCount += other.PartitionCount;
        SliceCount += other.SliceCount;
        SliceLength += other.SliceLength;

        for(int i = 0; i < SupportCounts.length; ++i)
        {
            SupportCounts[i] += other.SupportCounts[i];
        }
    }

    public String toString()
    {
        return format("partiions(%s) slices(%d totalLen=%d) reads(%d noVar=%s) readType(none=%d no=%d support=%d)",
                PartitionCount, SliceCount, SliceLength, ReadCount, NoVariantReadCount,
                SupportCounts[ReadMatchType.UNRELATED.ordinal()], SupportCounts[ReadMatchType.NO_SUPPORT.ordinal()],
                SupportCounts[ReadMatchType.SUPPORT.ordinal()]);
    }
}
