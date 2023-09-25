package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.Util.isNotProperReadPair;
import static com.hartwig.hmftools.sieve.annotate.Util.isSoftClipped;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class AnnotateStatistics
{
    public static final String CSV_HEADER =
            "PrimaryReadCount,PrimarySoftClippedCount,PrimaryImproperPairCount,PrimarySoftClippedANDImproperPairCount,SupplementaryCount";
    public static final String EMPTY_CSV_FRAGMENT = "NA,NA,NA,NA,NA";

    private long mPrimaryReadCount;
    private long mPrimarySoftClippedCount;
    private long mPrimaryImproperPairCount;
    private long mPrimarySoftClippedAndImproperPairCount;
    private long mSupplementaryCount;

    public AnnotateStatistics()
    {
        mPrimaryReadCount = 0;
        mPrimarySoftClippedCount = 0;
        mPrimaryImproperPairCount = 0;
        mPrimarySoftClippedAndImproperPairCount = 0;
        mSupplementaryCount = 0;
    }

    public void matchedRead(@NotNull final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            return;
        }

        mPrimaryReadCount++;

        boolean softClipped = isSoftClipped(read);
        boolean improperPair = isNotProperReadPair(read);

        if(softClipped && improperPair)
        {
            mPrimarySoftClippedAndImproperPairCount++;
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

    public String getCSVFragment()
    {
        return String.valueOf(mPrimaryReadCount) + ',' + mPrimarySoftClippedCount + ',' + mPrimaryImproperPairCount + ','
                + mPrimarySoftClippedAndImproperPairCount + ',' + mSupplementaryCount;
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

    public long getPrimarySoftClippedAndImproperPairCount()
    {
        return mPrimarySoftClippedAndImproperPairCount;
    }

    public long getSupplementaryCount()
    {
        return mSupplementaryCount;
    }
}
