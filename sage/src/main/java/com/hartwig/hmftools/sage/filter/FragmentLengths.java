package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSizeAbs;
import static com.hartwig.hmftools.sage.SageConstants.CHIMERIC_FRAGMENT_LENGTH_MAX;

import htsjdk.samtools.SAMRecord;

public class FragmentLengths
{
    private int mMaxAltLength;
    private int mNonAltCount;
    private long mNonAltLengthTotal;

    public FragmentLengths()
    {
        mMaxAltLength = 0;
        mNonAltCount = 0;
        mNonAltLengthTotal = 0;
    }

    public void processRead(final SAMRecord read, boolean supportsAlt)
    {
        int insertSize = inferredInsertSizeAbs(read);

        if(insertSize > CHIMERIC_FRAGMENT_LENGTH_MAX)
            return;

        if(supportsAlt)
        {
            mMaxAltLength = max(mMaxAltLength, insertSize);
        }
        else
        {
            ++mNonAltCount;
            mNonAltLengthTotal += insertSize;
        }
    }

    public int maxAltLength() { return mMaxAltLength; }

    public int nonAltCount() { return mNonAltCount; }
    public double averageNonAltLength() { return mNonAltCount > 0 ? mNonAltLengthTotal / (double)mNonAltCount : 0; }
}
