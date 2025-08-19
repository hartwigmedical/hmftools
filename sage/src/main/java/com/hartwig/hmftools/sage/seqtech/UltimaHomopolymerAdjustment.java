package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.calcTpBaseQual;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMRecord;

@VisibleForTesting
public class UltimaHomopolymerAdjustment extends UltimaQualModel
{
    private final int mHpStartIndex;
    private final int mHpEndIndex;
    private final int mRefAdjustCount; // bases to adjust to get back to the ref

    public UltimaHomopolymerAdjustment(final int hpStartIndex, final int hpEndIndex, final int refAdjustCount)
    {
        super(UltimaModelType.HOMOPOLYMER_ADJUSTMENT);

        mHpStartIndex = hpStartIndex;
        mHpEndIndex = hpEndIndex;
        mRefAdjustCount = refAdjustCount;
    }

    public byte calculateQual(final SAMRecord record, int varReadIndex)
    {
        return calcTpBaseQual(
                record, varReadIndex + mHpStartIndex, varReadIndex + mHpEndIndex, mRefAdjustCount);
    }

    public int hpStartIndex()
    {
        return mHpStartIndex;
    }
    public int hpEndIndex()
    {
        return mHpEndIndex;
    }
    public int refAdjustCount()
    {
        return mRefAdjustCount;
    }

    public String toString()
    {
        return format("adjust: hpIndex(%d-%d) adjust(%d)", mHpStartIndex, mHpEndIndex, mRefAdjustCount);
    }
}
