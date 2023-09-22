package com.hartwig.hmftools.sieve.annotate;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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

    // TODO(m_cooper): Use the original implementation.
    private static boolean isNotProperReadPair(@NotNull final SAMRecord read)
    {
        if(read.getReadUnmappedFlag())
        {
            return false;
        }

        // or a fragment length outside the observed distribution
        // TODO(m_cooper): Make this configurable.
        if(abs(read.getInferredInsertSize()) > 1000)
        {
            return true;
        }

        // an unmapped mate
        if(mateUnmapped(read))
        {
            return true;
        }

        if(read.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!read.getReferenceName().equals(read.getMateReferenceName()))
            {
                return true;
            }

            // inversion
            return read.getReadNegativeStrandFlag() == mateNegativeStrand(read);
        }

        return false;
    }

    public void matchedRead(@NotNull final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            mSupplementaryCount++;
            return;
        }

        mPrimaryReadCount++;

        final Cigar cigar = read.getCigar();
        final var it = cigar.getCigarElements().iterator();
        while(it.hasNext())
        {
            final CigarElement el = it.next();
            final CigarOperator op = el.getOperator();
            if(op == CigarOperator.SOFT_CLIP)
            {
                mPrimarySoftClippedCount++;
                break;
            }
        }

        if(isNotProperReadPair(read))
        {
            mPrimaryImproperPairCount++;
        }
    }

    public String getCSVFragment()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(mPrimaryReadCount);
        sb.append(',');
        sb.append(mPrimarySoftClippedCount);
        sb.append(',');
        sb.append(mPrimaryImproperPairCount);
        sb.append(',');
        sb.append(mSupplementaryCount);
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

    public long getSupplementaryCount()
    {
        return mSupplementaryCount;
    }

    public long getPrimaryImproperPairCount()
    {
        return mPrimaryImproperPairCount;
    }
}
