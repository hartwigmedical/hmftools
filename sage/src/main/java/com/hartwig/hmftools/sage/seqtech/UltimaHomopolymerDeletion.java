package com.hartwig.hmftools.sage.seqtech;

import static java.lang.String.format;

import static com.google.common.primitives.UnsignedBytes.max;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
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
    // the index of the bases relative to the variant that surround the deleted HP
    private final int mStraddleIndexStart;
    private final int mStraddleIndexEnd;
    private final byte mStraddleBaseStart;
    private final byte mStraddleBaseEnd;
    private final boolean mInCyclePosStrand;
    private final boolean mInCycleNegStrand;

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
            mStraddleBaseEnd = straddleBaseEnd; // the first ref base after the deleted bases

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
            byte altBase = (byte) variant.alt().charAt(0);

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
                return BQR_CACHE.maxRawQual();
        }
        else
        {
            if(!mInCyclePosStrand)
                return BQR_CACHE.maxRawQual();
        }

        final byte[] t0Values = extractT0Values(record);
        byte qual1 = safeQualLookup(t0Values, varReadIndex + mStraddleIndexStart);
        if(qual1 == ULTIMA_INVALID_QUAL)
        {
            return ULTIMA_INVALID_QUAL;
        }

        byte qual2 = safeQualLookup(t0Values, varReadIndex + mStraddleIndexEnd);
        if(qual2 == ULTIMA_INVALID_QUAL)
        {
            return ULTIMA_INVALID_QUAL;
        }

        return max(qual1, qual2);
    }

    public String toString()
    {
        return format("del: straddle(%d=%d %d=%d) cycle(pos=%s & neg=%s)",
                mStraddleIndexStart, mStraddleBaseStart, mStraddleIndexEnd, mStraddleBaseEnd,
                mInCyclePosStrand, mInCycleNegStrand);
    }
}
