package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.Util.isNotProperReadPair;
import static com.hartwig.hmftools.sieve.annotate.Util.isSoftClipped;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class AnnotateStatistics
{
    public static final String CSV_HEADER = "PrimaryReadCount,PrimarySoftClippedCount,PrimaryImproperPairCount,SupplementaryCount";
    public static final String EMPTY_CSV_FRAGMENT = "NA,NA,NA,NA";

    private long mPrimaryReadCount;
    private long mPrimarySoftClippedCount;
    private long mSupplementaryCount;
    private long mPrimaryImproperPairCount;

    public AnnotateStatistics()
    {
        mPrimaryReadCount = 0;
        mPrimarySoftClippedCount = 0;
        mSupplementaryCount = 0;
        mPrimaryImproperPairCount = 0;
    }

    public void matchedRead(@NotNull final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            return;
        }

        mPrimaryReadCount++;

        if(isSoftClipped(read))
        {
            mPrimarySoftClippedCount++;
        }

        if(isNotProperReadPair(read))
        {
            mPrimaryImproperPairCount++;
        }
    }

    public String getCSVFragment()
    {
        return String.valueOf(mPrimaryReadCount) + ',' + mPrimarySoftClippedCount + ',' + mPrimaryImproperPairCount + ','
                + mSupplementaryCount;
    }

    public long getPrimaryReadCount()
    {
        return mPrimaryReadCount;
    }

    public long getPrimarySoftClippedCount()
    {
        return mPrimarySoftClippedCount;
    }

    public long getSupplementaryCount()
    {
        return mSupplementaryCount;
    }

    public long getPrimaryImproperPairCount()
    {
        return mPrimaryImproperPairCount;
    }
}
