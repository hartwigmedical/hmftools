package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.google.common.primitives.UnsignedBytes.max;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.calcTpBaseQual;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isBaseInCycle;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.safeQualLookup;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class UltimaQualCalculator
{
    private static final int MAX_HOMOPOLYMER = 15; // temporarily since we don't fall back to MSI model properly
    private final RefGenomeInterface mRefGenome;

    public UltimaQualCalculator(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public UltimaQualModel buildContext(final SimpleVariant variant, final byte[] coreBases) // pass in MH or repeat info
    {
        int maxHomopolymerLength = Math.max(variant.ref().length(), MAX_HOMOPOLYMER);
        int refBaseEnd = variant.Position + maxHomopolymerLength + 1;

        // extract sufficient ref bases to set the context for most scenarios (only not for homopolymer transition)
        final byte[] refBases = mRefGenome.getBases(variant.Chromosome, variant.Position - 1, refBaseEnd);
        int refVarIndex = 1;

        if(variant.isIndel())
        {
            if(variant.isDelete() && isHomopolymerDeletion(variant, refBases))
            {
                return new HomopolymerDeletion(variant, refBases[2], coreBases[1], coreBases[2]);
            }

            int homopolymerTransitionIndex = findHomopolymerTransitionCandidate(variant);

            if(homopolymerTransitionIndex > 0)
            {
                int lowerHpEnd = variant.Position + homopolymerTransitionIndex - 1;
                int upperHpStart = lowerHpEnd + 1;
                int lowerRefBaseStart = lowerHpEnd - MAX_HOMOPOLYMER;
                int upperRefBaseEnd = upperHpStart + MAX_HOMOPOLYMER;
                final byte[] lowerRefBases = mRefGenome.getBases(variant.Chromosome, lowerRefBaseStart, lowerHpEnd);
                final byte[] upperRefBases = mRefGenome.getBases(variant.Chromosome, upperHpStart, upperRefBaseEnd);

                byte lowerHpBase = lowerRefBases[lowerRefBases.length - 1];
                int lowerHpLength = findHomopolymerLength(lowerRefBases, lowerHpBase, lowerRefBases.length - 1, false);
                byte upperHpBase = upperRefBases[0];
                int upperHpLength = findHomopolymerLength(upperRefBases, upperHpBase, 0, true);

                if(lowerHpLength > 2 && upperHpLength > 2 && variant.ref().length() + variant.position() - 1 <= upperRefBaseEnd)
                {
                    return new HomopolymerTransitionDeletion(variant, homopolymerTransitionIndex, lowerHpLength, upperHpLength);
                }
            }

            int homopolymerLength;
            if(variant.isInsert())
            {
                byte insertBase = (byte)variant.alt().charAt(1);
                homopolymerLength = findHomopolymerLength(refBases, insertBase, refVarIndex + 1, true);
            }
            else
            {
                // the HP search start 1 base after the variant's ref position
                byte hpBase = refBases[refVarIndex + 1];
                homopolymerLength = findHomopolymerLength(refBases, hpBase, refVarIndex + 1, true);
            }

            boolean singleRefBase = variant.ref().substring(1).chars().distinct().count() <= 1;
            boolean singleAltBase = variant.alt().substring(1).chars().distinct().count() <= 1;
            if(homopolymerLength <= MAX_HOMOPOLYMER && singleRefBase && singleAltBase)
            {
                int homopolymerStartIndex = 1;
                // adjust the HP length by the diff observed in the variant
                int homopolymerEndIndex = homopolymerStartIndex + homopolymerLength - 1 + variant.indelLength();
                int refAdjustCount = -variant.indelLength();

                return new HomopolymerAdjustment(homopolymerStartIndex, homopolymerEndIndex, refAdjustCount);
            }
            else
            {
                return new MicrosatelliteAdjustment();
            }
        }
        else if(variant.isSNV())
        {
            return new SnvMnv(variant, refBases, refVarIndex, coreBases[0], coreBases[2]);
        }

        return null;
    }

    private static int findHomopolymerTransitionCandidate(final SimpleVariant variant)
    {
        // returns the index within the deleted ref bases of a transition, otherwise 0
        if(!variant.isDelete() || variant.ref().length() < 3)
            return 0;

        char firstDelBase = variant.ref().charAt(1);

        for(int i = 2; i < variant.ref().length(); ++i)
        {
            if(variant.ref().charAt(i) != firstDelBase)
                return i;
        }

        return 0;
    }

    private static boolean isHomopolymerDeletion(final SimpleVariant variant, final byte[] refBases)
    {
        byte refBase = refBases[1];
        byte deletedBase = refBases[2];
        int delLength = abs(variant.indelLength());
        byte postDelRefBase = refBases[1 + variant.Ref.length()];

        if(refBase == deletedBase || postDelRefBase == deletedBase)
            return false;

        // for a DEL of 2 bases this will check the 3rd base, being then second deleted base, and then stop
        for(int i = 3; i <= 1 + delLength; ++i)
        {
            if(refBases[i] != deletedBase)
                return false;
        }

        return true;
    }

    private static int findHomopolymerLength(final byte[] refBases, final byte compareBase, int startIndex, boolean searchUp)
    {
        // byte repeatBase = refBases[startIndex];
        int repeatCount = 0;

        // int i = startIndex + (searchUp ? 1 : -1);
        int i = startIndex;

        while(i >= 0 && i < refBases.length)
        {
            if(refBases[i] != compareBase)
                return repeatCount;

            ++repeatCount;

            if(searchUp)
                ++i;
            else
                --i;
        }

        return repeatCount;
    }

    private class HomopolymerAdjustment extends UltimaQualModel
    {
        private final int mHpStartIndex;
        private final int mHpEndIndex;
        private final int mRefAdjustCount;

        public HomopolymerAdjustment(final int hpStartIndex, final int hpEndIndex, final int refAdjustCount)
        {
            super(UltimaModelType.HOMOPOLYMER_ADJUSTMENT);

            // TODO: Do we keep this around after testing?
//            if(hpEndIndex < hpStartIndex)
//            {
//                throw new IllegalArgumentException("Homopolymer end index is before start index.");
//            }

            mHpStartIndex = hpStartIndex;
            mHpEndIndex = hpEndIndex;
            mRefAdjustCount = refAdjustCount;
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex)
        {
            return calcTpBaseQual(
                    record, varReadIndex + mHpStartIndex, varReadIndex + mHpEndIndex, mRefAdjustCount);
        }

        public String toString()
        {
            return format("adjust: hpIndex(%d-%d) adjust(%d)", mHpStartIndex, mHpEndIndex, mRefAdjustCount);
        }
    }

    private class HomopolymerDeletion extends UltimaQualModel
    {
        // the index of the bases relative to the variant that surround the deleted HP
        private final int mStraddleIndexStart;
        private final int mStraddleIndexEnd;
        private final byte mStraddleBaseStart;
        private final byte mStraddleBaseEnd;
        private final boolean mInCyclePosStrand;
        private final boolean mInCycleNegStrand;

        public HomopolymerDeletion(final SimpleVariant variant, final byte deletedHpBase, final byte straddleBaseStart, final byte straddleBaseEnd)
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
            }
            else
            {
                // used for SNVs where the variant's ref base is considered deleted
                mStraddleIndexStart = 0;
                mStraddleIndexEnd = 0;

                // the 2 bases surrounding the SNV
                mStraddleBaseStart = straddleBaseStart;
                mStraddleBaseEnd = straddleBaseEnd;

                deletedBase = deletedHpBase;
            }

            mInCyclePosStrand = isBaseInCycle(mStraddleBaseStart, mStraddleBaseEnd, deletedBase);

            byte revDeletedBase = swapDnaBase(deletedBase);
            mInCycleNegStrand = isBaseInCycle(swapDnaBase(mStraddleBaseEnd), swapDnaBase(mStraddleBaseStart), revDeletedBase);
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex)
        {
            if(record.getReadNegativeStrandFlag())
            {
                if(!mInCycleNegStrand)
                    return ULTIMA_MAX_QUAL;
            }
            else
            {
                if(!mInCyclePosStrand)
                    return ULTIMA_MAX_QUAL;
            }

            final byte[] t0Values = record.getStringAttribute(T0_TAG).getBytes();
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

            return (byte) (max(qual1, qual2) - 33);
        }

        public String toString()
        {
            return format("del: straddle(%d=%d %d=%d) cycle(pos=%s & neg=%s)",
                    mStraddleIndexStart, mStraddleBaseStart, mStraddleIndexEnd, mStraddleBaseEnd,
                    mInCyclePosStrand, mInCycleNegStrand);
        }
    }

    private class HomopolymerTransitionDeletion extends UltimaQualModel
    {
        private final int mLowerHpStartIndex;
        private final int mLowerHpEndIndex;
        private final int mLowerRefAdjustCount;
        private final int mUpperHpStartIndex;
        private final int mUpperHpEndIndex;
        private final int mUpperRefAdjustCount;

        public HomopolymerTransitionDeletion(
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

            return (byte)(lowerQual + upperQual);
        }

        public String toString()
        {
            return format("transition: lower(%d-%d %d) upper(%d-%d %d)",
                    mLowerHpStartIndex, mLowerHpEndIndex, mLowerRefAdjustCount,
                    mUpperHpStartIndex, mUpperHpEndIndex, mLowerRefAdjustCount);
        }
    }

    private class MicrosatelliteAdjustment extends UltimaQualModel
    {
        public MicrosatelliteAdjustment()
        {
            super(UltimaModelType.MICROSAT_ADJUSTMENT);
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return UltimaBamUtils.ULTIMA_MAX_QUAL; }
    }

    private class SnvMnv extends UltimaQualModel
    {
        private final HomopolymerAdjustment mLeftAdjust;
        private final HomopolymerAdjustment mRightAdjust;
        private final HomopolymerDeletion mLeftDeletion;
        private final HomopolymerDeletion mRightDeletion;

        public SnvMnv(
                final SimpleVariant variant, final byte[] refBases, final int refVarIndex, final byte leftReadBase, final byte rightReadBase)
        {
            super(UltimaModelType.SNV);

            byte refBase = (byte)variant.ref().charAt(0);
            byte altBase = (byte)variant.alt().charAt(0);
            byte leftBase = refBases[refVarIndex - 1];
            byte rightBase = refBases[refVarIndex + 1];

            // SNVs and MNVs
            if(leftBase == rightBase && (leftBase == refBase || leftBase == altBase))
            {
                // both straddling bases match either the ref or alt, so no need to check for HP adjustment
                mLeftAdjust = null;
                mRightAdjust = null;
                mLeftDeletion = null;
                mRightDeletion = null;
                return;
            }

            // one side will be a delete (either partial or complete), the other a 1-base insert
            HomopolymerAdjustment leftAdjust = null;
            HomopolymerAdjustment rightAdjust = null;
            HomopolymerDeletion leftDeletion = null;
            HomopolymerDeletion rightDeletion = null;

            boolean leftMatchesRef = refBases[refVarIndex - 1] == refBases[refVarIndex];
            boolean leftMatchesAlt = refBases[refVarIndex - 1] == altBase;
            boolean rightMatchesRef = refBases[refVarIndex + 1] == refBases[refVarIndex];
            boolean rightMatchesAlt = refBases[refVarIndex + 1] == altBase;
            boolean leftSideDelete = (leftMatchesRef && !leftMatchesAlt) || rightMatchesAlt;

            int lowerHpLength = 0;
            if(leftMatchesAlt || leftMatchesRef)
            {
                int lowerRefBaseEnd = leftMatchesRef ? variant.Position : variant.Position - 1;
                int lowerRefBaseStart = lowerRefBaseEnd - MAX_HOMOPOLYMER;
                final byte[] lowerRefBases = mRefGenome.getBases(variant.Chromosome, lowerRefBaseStart, lowerRefBaseEnd);
                byte hpBase = lowerRefBases[lowerRefBases.length - 1];
                lowerHpLength = findHomopolymerLength(lowerRefBases, hpBase, lowerRefBases.length - 1, false);
            }

            if(leftMatchesAlt || !leftSideDelete)
            {
                // is a HP expansion of 1 base on the left - find DP length start with the left position
                int newHpLength = lowerHpLength + 1;

                int leftHpEndIndex = 0;
                int leftHpStartIndex = leftHpEndIndex - newHpLength + 1;
                leftAdjust = new HomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, -1);
            }
            else if(leftMatchesRef)
            {
                // is a HP contraction of 1 base on the left, starting with the variant's position
                int newLowerHpLength = lowerHpLength - 1; // reflecting the new length of the HP after the delete

                int leftHpEndIndex = -1;
                int leftHpStartIndex = leftHpEndIndex - newLowerHpLength + 1;
                leftAdjust = new HomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, 1);
            }
            else
            {
                // HP deletion
                leftDeletion = new HomopolymerDeletion(variant, refBases[1], leftReadBase, rightReadBase);
            }

            int upperHpLength = 0;
            if(rightMatchesAlt || rightMatchesRef)
            {
                int compareStartIndex = rightMatchesRef ? 1 : 2; // 1 being the variant's base, 2 being the one after
                byte hpBase = refBases[2];
                upperHpLength = findHomopolymerLength(refBases, hpBase, compareStartIndex, true);
            }

            if(rightMatchesAlt || leftSideDelete)
            {
                // is a HP expansion of 1 base on the right
                int newHpLength = upperHpLength + 1;

                int rightHpStartIndex = 0;
                int rightHpEndIndex = rightHpStartIndex + newHpLength - 1;
                rightAdjust = new HomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, -1);
            }
            else if(rightMatchesRef)
            {
                // is a HP contraction of 1 base on the right, starting with the variant's position
                int newHpLength = upperHpLength - 1;

                int rightHpStartIndex = 1;
                int rightHpEndIndex = rightHpStartIndex + newHpLength - 1;
                rightAdjust = new HomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, 1);
            }
            else
            {
                rightDeletion = new HomopolymerDeletion(variant, refBases[1], leftReadBase, rightReadBase);
            }

            mLeftAdjust = leftAdjust;
            mRightAdjust = rightAdjust;
            mLeftDeletion = leftDeletion;
            mRightDeletion = rightDeletion;
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex)
        {
            if(mLeftAdjust == null && mLeftDeletion == null)
                return ULTIMA_MAX_QUAL;

            if(mLeftDeletion != null || mRightDeletion != null)
            {
                // any HP deletion can just use the SNV base itself
                final byte[] t0Values = record.getStringAttribute(T0_TAG).getBytes();
                return safeQualLookup(t0Values, varReadIndex);
            }

            int leftQual = mLeftAdjust != null ?
                    mLeftAdjust.calculateQual(record, varReadIndex) : mLeftDeletion.calculateQual(record, varReadIndex);

            int rightQual = mRightAdjust != null ?
                    mRightAdjust.calculateQual(record, varReadIndex) : mRightDeletion.calculateQual(record, varReadIndex);

            int combinedQual = min(leftQual + rightQual, ULTIMA_MAX_QUAL);

            return (byte)combinedQual;
        }

        public String toString()
        {
            if(mLeftDeletion == null && mLeftAdjust == null)
                return "non-HP";

            return format("left(%s) right(%s)",
                    mLeftAdjust != null ? mLeftAdjust : mLeftDeletion, mRightAdjust != null ? mRightAdjust : mRightDeletion);
        }
    }
}
