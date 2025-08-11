package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.calcTpBaseQual;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMRecord;

@VisibleForTesting
public class UltimaHomopolymerAdjustment extends UltimaQualModel
{
    private final int mHpStartIndex;
    private final int mHpEndIndex;
    private final int mRefAdjustCount;

    public UltimaHomopolymerAdjustment(final int hpStartIndex, final int hpEndIndex, final int refAdjustCount)
    {
        super(UltimaModelType.HOMOPOLYMER_ADJUSTMENT);

        // TODO: Do we keep this around after testing?
        if(hpEndIndex < hpStartIndex)
        {
            throw new IllegalArgumentException("Homopolymer end index is before start index.");
        }

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
