package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.BQR_CACHE;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.calcTpBaseQual;

import com.hartwig.hmftools.common.variant.SimpleVariant;

import htsjdk.samtools.SAMRecord;

class UltimaHomopolymerTransitionDeletion extends UltimaQualModel
{
    private final int mLowerHpStartIndex;
    private final int mLowerHpEndIndex;
    private final int mLowerRefAdjustCount;
    private final int mUpperHpStartIndex;
    private final int mUpperHpEndIndex;
    private final int mUpperRefAdjustCount;

    public UltimaHomopolymerTransitionDeletion(
            final SimpleVariant variant, final int homopolymerTransitionIndex, final int lowerHpLength, final int upperHpLength)
    {
        super(UltimaModelType.HOMOPOLYMER_TRANSITION);

        int delLength = abs(variant.indelLength());
        int lowerDelLength = homopolymerTransitionIndex - 1;
        int upperDelLength = delLength - lowerDelLength;

        // the ref base always matches the lower HP base, and the ref + 1 base the upper HP base
        mLowerHpEndIndex = 0;
        mLowerHpStartIndex = mLowerHpEndIndex - lowerHpLength + 1 + lowerDelLength;
        mLowerRefAdjustCount = lowerDelLength;

        mUpperHpStartIndex = mLowerHpEndIndex + 1;
        mUpperHpEndIndex = mUpperHpStartIndex + upperHpLength - 1 - upperDelLength;
        mUpperRefAdjustCount = upperDelLength;
    }

    public byte calculateQual(final SAMRecord record, int varReadIndex)
    {
        byte lowerQual = calcTpBaseQual(
                record, varReadIndex + mLowerHpStartIndex, varReadIndex + mLowerHpEndIndex, mLowerRefAdjustCount);

        if(lowerQual == ULTIMA_INVALID_QUAL)
            return ULTIMA_INVALID_QUAL;

        byte upperQual = calcTpBaseQual(
                record, varReadIndex + mUpperHpStartIndex, varReadIndex + mUpperHpEndIndex, mUpperRefAdjustCount);

        if(upperQual == ULTIMA_INVALID_QUAL)
            return ULTIMA_INVALID_QUAL;

        // cap at higher of BQR quals
        int lowerHpLength = mLowerHpEndIndex - mLowerHpStartIndex + 1;
        char lowerHpBase = (char)record.getReadBases()[varReadIndex + mLowerHpStartIndex];
        int lowerBqrQual = BQR_CACHE.getTpRecalibratedQual(lowerHpLength, lowerHpBase, false);

        int upperHpLength = mUpperHpEndIndex - mUpperHpStartIndex + 1;
        char upperHpBase = (char)record.getReadBases()[varReadIndex + mUpperHpStartIndex];
        int upperBqrQual = BQR_CACHE.getTpRecalibratedQual(upperHpLength, upperHpBase, false);

        int maxBqrQual = max(lowerBqrQual, upperBqrQual);

        return (byte)min(lowerQual + upperQual, maxBqrQual);
    }

    public String toString()
    {
        return format("transition: lower(%d-%d %d) upper(%d-%d %d)",
                mLowerHpStartIndex, mLowerHpEndIndex, mLowerRefAdjustCount,
                mUpperHpStartIndex, mUpperHpEndIndex, mLowerRefAdjustCount);
    }
}
