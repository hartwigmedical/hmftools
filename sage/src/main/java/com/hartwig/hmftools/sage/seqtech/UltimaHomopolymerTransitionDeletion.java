package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_TP;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_0_BOOST;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.calcTpBaseQual;

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
        {
            return ULTIMA_INVALID_QUAL;
        }

        byte upperQual = calcTpBaseQual(
                record, varReadIndex + mUpperHpStartIndex, varReadIndex + mUpperHpEndIndex, mUpperRefAdjustCount);

        if(upperQual == ULTIMA_INVALID_QUAL)
        {
            return ULTIMA_INVALID_QUAL;
        }

        if(lowerQual == ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST || upperQual == ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST)
        {
            return ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST;
        }

        return (byte) min(lowerQual + upperQual, ULTIMA_MAX_QUAL_TP);
    }

    public String toString()
    {
        return format("transition: lower(%d-%d %d) upper(%d-%d %d)",
                mLowerHpStartIndex, mLowerHpEndIndex, mLowerRefAdjustCount,
                mUpperHpStartIndex, mUpperHpEndIndex, mLowerRefAdjustCount);
    }
}
