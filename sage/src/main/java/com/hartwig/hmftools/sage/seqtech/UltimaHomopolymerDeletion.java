package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.google.common.primitives.UnsignedBytes.max;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isBaseInCycle;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.BQR_CACHE;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.safeQualLookup;

import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import htsjdk.samtools.SAMRecord;

class UltimaHomopolymerDeletion extends UltimaQualModel
{
    // the index and read base straddling the variant's deleted HP
    private final int mStraddleIndexStart;
    private final int mStraddleIndexEnd;

    private final byte mStraddleBaseStart;
    private final byte mStraddleBaseEnd;

    // whether the delete is in-cycle when read from each direction
    private final boolean mInCyclePosStrand;
    private final boolean mInCycleNegStrand;

    public UltimaHomopolymerDeletion(
            final int straddleIndexStart, final int straddleIndexEnd, final byte straddleBaseStart,
            final byte straddleBaseEnd, final boolean inCyclePosStrand, final boolean inCycleNegStrand)
    {
        super(UltimaModelType.HOMOPOLYMER_DELETION);

        mStraddleIndexStart = straddleIndexStart;
        mStraddleIndexEnd = straddleIndexEnd;
        mStraddleBaseStart = straddleBaseStart;
        mStraddleBaseEnd = straddleBaseEnd;
        mInCyclePosStrand = inCyclePosStrand;
        mInCycleNegStrand = inCycleNegStrand;
    }

    public static UltimaHomopolymerDeletion fromDelete(final byte deletedBase, final byte refBase, final byte postDeleteRefBase)
    {
        boolean inCyclePosStrand = isBaseInCycle(refBase, postDeleteRefBase, deletedBase);
        boolean inCycleNegStrand = isBaseInCycle(swapDnaBase(postDeleteRefBase), swapDnaBase(refBase), swapDnaBase(deletedBase));

        return new UltimaHomopolymerDeletion(
                0, 1, refBase, postDeleteRefBase, inCyclePosStrand, inCycleNegStrand);
    }

    public static UltimaHomopolymerDeletion fromMnv(
            final byte deletedRefBase, final byte newAltBase, final byte readBaseBefore, final byte readBaseAfter)
    {
        // for an SNV, the deleted ref will be the ref of the SNV, and similarly the new alt will be the alt
        // otherwise it will be the same for a specific position within an MNV
        int posCycleCountRef = UltimaBamUtils.cycleCount(readBaseBefore, readBaseAfter, deletedRefBase);
        int posCycleCountAlt = UltimaBamUtils.cycleCount(readBaseBefore, readBaseAfter, newAltBase);
        boolean inCyclePosStrand = posCycleCountRef == posCycleCountAlt;

        int negCycleCountRef = UltimaBamUtils.cycleCount(swapDnaBase(readBaseAfter), swapDnaBase(readBaseBefore), swapDnaBase(deletedRefBase));
        int negCycleCountAlt = UltimaBamUtils.cycleCount(swapDnaBase(readBaseAfter), swapDnaBase(readBaseBefore), swapDnaBase(newAltBase));
        boolean inCycleNegStrand = negCycleCountRef == negCycleCountAlt;

        return new UltimaHomopolymerDeletion(
                0, 0, readBaseBefore, readBaseAfter, inCyclePosStrand, inCycleNegStrand);
    }

    @Deprecated
    public UltimaHomopolymerDeletion(
            final SimpleVariant variant, final byte deletedHpBase, final byte straddleBaseStart, final byte straddleBaseEnd)
    {
        super(UltimaModelType.HOMOPOLYMER_DELETION);

        byte deletedBase;

        if(variant.isDelete())
        {
            mStraddleIndexStart = 0;
            mStraddleIndexEnd = 1;
            mStraddleBaseStart = straddleBaseStart; // the variant position ref base
            mStraddleBaseEnd = straddleBaseEnd; // the first read base after the deleted bases

            deletedBase = deletedHpBase;

            mInCyclePosStrand = isBaseInCycle(mStraddleBaseStart, mStraddleBaseEnd, deletedBase);

            byte revDeletedBase = swapDnaBase(deletedBase);
            mInCycleNegStrand = isBaseInCycle(swapDnaBase(mStraddleBaseEnd), swapDnaBase(mStraddleBaseStart), revDeletedBase);
        }
        else  // used for SNVs where the variant's ref base is considered deleted
        {
            // any HP deletion can just use the SNV base itself
            mStraddleIndexStart = 0;
            mStraddleIndexEnd = 0;

            // the 2 bases surrounding the SNV
            mStraddleBaseStart = straddleBaseStart;
            mStraddleBaseEnd = straddleBaseEnd;

            deletedBase = deletedHpBase;
            byte altBase = (byte)variant.alt().charAt(0);

            int posCycleCountRef = UltimaBamUtils.cycleCount(straddleBaseStart, mStraddleBaseEnd, deletedBase);
            int posCycleCountAlt = UltimaBamUtils.cycleCount(straddleBaseStart, mStraddleBaseEnd, altBase);
            mInCyclePosStrand = posCycleCountRef == posCycleCountAlt;

            int negCycleCountRef = UltimaBamUtils.cycleCount(
                    swapDnaBase(mStraddleBaseEnd), swapDnaBase(mStraddleBaseStart), swapDnaBase(deletedHpBase));

            int negCycleCountAlt = UltimaBamUtils.cycleCount(
                    swapDnaBase(mStraddleBaseEnd), swapDnaBase(mStraddleBaseStart), swapDnaBase(altBase));

            mInCycleNegStrand = negCycleCountRef == negCycleCountAlt;
        }
    }

    public byte calculateQual(final SAMRecord record, int varReadIndex)
    {
        if(record.getReadNegativeStrandFlag())
        {
            if(!mInCycleNegStrand)
                return BQR_CACHE.outOfCycleT0Qual();
        }
        else
        {
            if(!mInCyclePosStrand)
                return BQR_CACHE.outOfCycleT0Qual();
        }

        final byte[] t0Values = extractT0Values(record);
        byte qual1 = safeQualLookup(t0Values, varReadIndex + mStraddleIndexStart);

        qual1 = BQR_CACHE.calcT0RecalibratedQual(qual1, record, varReadIndex + mStraddleIndexStart);

        if(qual1 == ULTIMA_INVALID_QUAL)
            return ULTIMA_INVALID_QUAL;

        byte qual2;
        if(mStraddleIndexEnd == mStraddleIndexStart)
        {
            qual2 = qual1;
        }
        else
        {
            qual2 = safeQualLookup(t0Values, varReadIndex + mStraddleIndexEnd);
            qual2 = BQR_CACHE.calcT0RecalibratedQual(qual2, record, varReadIndex + mStraddleIndexEnd);

            if(qual2 == ULTIMA_INVALID_QUAL)
                return ULTIMA_INVALID_QUAL;
        }

        return maxQual(qual1, qual2);
    }

    public String toString()
    {
        return format("del: straddle(%d=%d %d=%d) cycle(pos=%s & neg=%s)",
                mStraddleIndexStart, mStraddleBaseStart, mStraddleIndexEnd, mStraddleBaseEnd,
                mInCyclePosStrand, mInCycleNegStrand);
    }
}
