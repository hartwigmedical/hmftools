package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.Util.isNotProperReadPair;
import static com.hartwig.hmftools.sieve.annotate.Util.isSoftClipped;

import htsjdk.samtools.SAMRecord;

public class AnnotateStatistics
{
    // Bucket MapQs into 0-9, 10-19, ..., 40-49, >= 50.
    private static final int MAPQ_BUCKET_COUNT = 6;
    public static final String TSV_HEADER = getTSVHeader();

    private final long[] mPrimaryMapQBuckets;
    private final long[] mSupplementaryMapQBuckets;

    private long mPrimaryReadCount;
    private long mPrimarySoftClippedCount;
    private long mPrimaryImproperPairCount;
    private long mPrimarySoftClippedOrImproperPairCount;
    private long mSupplementaryCount;

    public AnnotateStatistics()
    {
        mPrimaryReadCount = 0;
        mPrimarySoftClippedCount = 0;
        mPrimaryImproperPairCount = 0;
        mPrimarySoftClippedOrImproperPairCount = 0;
        mSupplementaryCount = 0;

        mPrimaryMapQBuckets = new long[MAPQ_BUCKET_COUNT];
        for(int i = 0; i < mPrimaryMapQBuckets.length; i++)
        {
            mPrimaryMapQBuckets[i] = 0;
        }

        mSupplementaryMapQBuckets = new long[MAPQ_BUCKET_COUNT];
        for(int i = 0; i < mSupplementaryMapQBuckets.length; i++)
        {
            mSupplementaryMapQBuckets[i] = 0;
        }
    }

    private static String getTSVHeader()
    {
        final StringBuilder sb = new StringBuilder();
        sb.append("PrimaryReadCount");
        sb.append('\t');
        sb.append("PrimarySoftClippedCount");
        sb.append('\t');
        sb.append("PrimaryImproperPairCount");
        sb.append('\t');
        sb.append("PrimarySoftClippedORImproperPairCount");
        sb.append('\t');
        sb.append("SupplementaryCount");

        for(int i = 0; i < MAPQ_BUCKET_COUNT; i++)
        {
            final int mapQStart = 10 * i;
            final int mapQEnd = mapQStart + 9;
            if(i < MAPQ_BUCKET_COUNT - 1)
            {
                sb.append("\tPrimary MAPQ " + mapQStart + '-' + mapQEnd);
            }
            else
            {
                sb.append("\tPrimary MAPQ >=" + mapQStart);
            }
        }

        for(int i = 0; i < MAPQ_BUCKET_COUNT; i++)
        {
            final int mapQStart = 10 * i;
            final int mapQEnd = mapQStart + 9;
            if(i < MAPQ_BUCKET_COUNT - 1)
            {
                sb.append("\tSupplementary MAPQ " + mapQStart + '-' + mapQEnd);
            }
            else
            {
                sb.append("\tSupplementary MAPQ >=" + mapQStart);
            }
        }

        return sb.toString();
    }

    private static int getMapQBucketIndex(int mapQ)
    {
        return Math.min(mapQ / 10, MAPQ_BUCKET_COUNT - 1);
    }

    public void matchedRead(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            mSupplementaryMapQBuckets[getMapQBucketIndex(read.getMappingQuality())]++;
            return;
        }

        mPrimaryReadCount++;
        mPrimaryMapQBuckets[getMapQBucketIndex(read.getMappingQuality())]++;

        boolean softClipped = isSoftClipped(read);
        boolean improperPair = isNotProperReadPair(read);

        if(softClipped || improperPair)
        {
            mPrimarySoftClippedOrImproperPairCount++;
        }

        if(softClipped)
        {
            mPrimarySoftClippedCount++;
        }

        if(improperPair)
        {
            mPrimaryImproperPairCount++;
        }
    }

    public String getTSVFragment()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(String.valueOf(mPrimaryReadCount) + '\t' + mPrimarySoftClippedCount + '\t' + mPrimaryImproperPairCount + '\t'
                + mPrimarySoftClippedOrImproperPairCount + '\t' + mSupplementaryCount);

        for(int i = 0; i < mPrimaryMapQBuckets.length; i++)
        {
            sb.append('\t');
            sb.append(mPrimaryMapQBuckets[i]);
        }

        for(int i = 0; i < mSupplementaryMapQBuckets.length; i++)
        {
            sb.append('\t');
            sb.append(mSupplementaryMapQBuckets[i]);
        }

        return sb.toString();
    }

    public long getPrimaryReadCount()
    {
        return mPrimaryReadCount;
    }

    public long getPrimarySoftClippedCount()
    {
        return mPrimarySoftClippedCount;
    }

    public long getPrimaryImproperPairCount()
    {
        return mPrimaryImproperPairCount;
    }

    public long getPrimarySoftClippedOrImproperPairCount()
    {
        return mPrimarySoftClippedOrImproperPairCount;
    }

    public long getSupplementaryCount()
    {
        return mSupplementaryCount;
    }
}
