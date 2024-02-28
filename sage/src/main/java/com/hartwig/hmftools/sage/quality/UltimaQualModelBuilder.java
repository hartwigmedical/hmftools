package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;

import static com.google.common.primitives.UnsignedBytes.max;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.calcTpBaseQual;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isBaseInCycle;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelBuilder
{

    private static final int MAX_HOMOPOLYMER = 8;

    private final RefGenomeInterface mRefGenome;

    public UltimaQualModelBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public UltimaQualModel buildContext(final SimpleVariant variant) // pass in MH or repeat info
    {
        int positionEnd = variant.Position + variant.ref().length();
        int refBaseEnd = variant.Position + MAX_HOMOPOLYMER + 1;

        // extract sufficient ref bases to set the context for each scenario
        final byte[] refBases = mRefGenome.getBases(variant.Chromosome, variant.Position - 1, refBaseEnd);

        if(variant.isIndel())
        {
            if(variant.isDelete() && isHomopolymerDeletion(variant, refBases))
            {
                return new HomopolymerDeletion(variant, refBases);
            }

            int homopolymerTransitionIndex = findHomopolymerTransitionCandidate(variant);

            if(homopolymerTransitionIndex > 0)
            {
                int lowerHpStart = variant.Position + homopolymerTransitionIndex;
                int upperHpStart = lowerHpStart + 1;
                int lowerRefBaseStart = lowerHpStart - MAX_HOMOPOLYMER;
                int upperRefBaseEnd = upperHpStart + MAX_HOMOPOLYMER;
                final byte[] lowerRefBases = mRefGenome.getBases(variant.Chromosome, lowerRefBaseStart, lowerHpStart);
                final byte[] upperRefBases = mRefGenome.getBases(variant.Chromosome, upperHpStart, upperRefBaseEnd);
                int lowerHomopolymerLength = findHomopolymerLength(lowerRefBases, lowerRefBases.length - 1, false);
                int upperHomopolymerLength = findHomopolymerLength(upperRefBases, 0, true);

                if(lowerHomopolymerLength > 2 && upperHomopolymerLength > 2)
                {
                    int delLength = abs(variant.indelLength());
                    int upperDelLength = delLength - homopolymerTransitionIndex;
                    int lowerDelLength = delLength - delLength;

                    int lowerHpIndexEnd = homopolymerTransitionIndex;
                    int lowerHpIndexStart = lowerHpIndexEnd - lowerHomopolymerLength + 1;
                    int upperHpIndexStart = lowerHpIndexEnd + 1;
                    int upperHpIndexEnd = upperHpIndexStart + upperHomopolymerLength - 1;

                    return new HomopolymerTransitionDeletion(
                            lowerHpIndexStart, lowerHpIndexEnd, -lowerDelLength, upperHpIndexEnd, -upperDelLength);
                }
            }

            int homopolymerLength = findHomopolymerLength(refBases);

            if(homopolymerLength <= MAX_HOMOPOLYMER)
            {
                int homopolymerStartIndex = 1;
                int homopolymerEndIndex = homopolymerStartIndex + homopolymerLength - 1;
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
            // SNVs and MNVs
            int homopolymerStartIndex = 0;
            int homopolymerEndIndex = 0;
            int refAdjustCount = 0;

            if(refBases[0] == refBases[2])
            {
                // ref and alt match, no need to check for HP adjustment
            }
            else if(refBases[0] == refBases[1] || refBases[1] == refBases[2])
            {
                // like a 1-base insertion of a base on the lower side with a 1-base deletion on the upper, or vice versa
                if(refBases[0] == refBases[1])
                {
                    // search down
                    int lowerRefBaseStart = variant.Position - MAX_HOMOPOLYMER;
                    final byte[] lowerRefBases = mRefGenome.getBases(variant.Chromosome, lowerRefBaseStart, variant.Position);
                    int lowerHomopolymerLength = findHomopolymerLength(lowerRefBases, lowerRefBases.length - 1, false);
                    homopolymerEndIndex = 0;
                    homopolymerStartIndex = homopolymerEndIndex - lowerHomopolymerLength + 1;
                }
                else
                {
                    int upperHomopolymerLength = findHomopolymerLength(refBases, 1, true);
                    homopolymerStartIndex = 0;
                    homopolymerEndIndex = homopolymerStartIndex + upperHomopolymerLength - 1;
                }
            }

            return new SnvMnv(variant, refBases, homopolymerStartIndex, homopolymerEndIndex, refAdjustCount);
        }
        else
        {
            return null;
        }
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
        int nextUndeletedBase = variant.Ref.length() + 2;
        return refBases[1] == refBases[nextUndeletedBase]; // will be checking bases 1 & 3 for a 1-base DEL
    }

    private static int findHomopolymerLength(final byte[] refBases) { return findHomopolymerLength(refBases, 2, true); }

    private static int findHomopolymerLength(final byte[] refBases, int startIndex, boolean searchUp)
    {
        byte repeatBase = refBases[startIndex]; // the base afer the variant's position
        int repeatCount = 1;

        int i = startIndex + (searchUp ? 1 : -1);

        while(i >= 0 && i < refBases.length)
        {
            if(refBases[i] != repeatBase)
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
            super(UltimaVariantModel.HOMOPOLYMER_ADJUSTMENT);
            mHpStartIndex = hpStartIndex;
            mHpEndIndex = hpEndIndex;
            mRefAdjustCount = refAdjustCount;
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex)
        {
            return calcTpBaseQual(
                    record, varReadIndex + mHpStartIndex, varReadIndex + mHpEndIndex, mRefAdjustCount);
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

        public HomopolymerDeletion(final SimpleVariant variant, final byte[] refBases)
        {
            super(UltimaVariantModel.HOMOPOLYMER_DELETION);
            mStraddleIndexStart = 0;
            mStraddleIndexEnd = variant.Ref.length();
            mStraddleBaseStart = refBases[1];
            mStraddleBaseEnd = refBases[1 + mStraddleIndexEnd];

            char deletedBase = variant.Ref.charAt(1);
            mInCyclePosStrand = isBaseInCycle(mStraddleBaseStart, mStraddleBaseEnd, (byte)deletedBase);

            char revDeletedBase = swapDnaBase(deletedBase);
            mInCycleNegStrand = isBaseInCycle(swapDnaBase(mStraddleBaseEnd), swapDnaBase(mStraddleBaseStart), (byte)revDeletedBase);
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

            final String t0Values = record.getStringAttribute(T0_TAG);
            byte qual1 = (byte)t0Values.charAt(mStraddleIndexStart);
            byte qual2 = (byte)t0Values.charAt(mStraddleIndexStart);

            return max(qual1, qual2);
        }
    }

    private class HomopolymerTransitionDeletion extends UltimaQualModel
    {
        private final int mLowerHomopolymerStartIndex;
        private final int mLowerHomopolymerEndIndex;
        private final int mLowerRefAdjustCount;
        private final int mUpperHomopolymerStartIndex;
        private final int mUpperHomopolymerEndIndex;
        private final int mUpperRefAdjustCount;

        public HomopolymerTransitionDeletion(
                final int lowerHpStartIndex, final int lowerHpEndIndex, final int lowerRefAdjustCount,
                final int upperHpEndIndex, final int upperRefAdjustCount)
        {
            super(UltimaVariantModel.HOMOPOLYMER_ADJUSTMENT);
            mLowerHomopolymerStartIndex = lowerHpStartIndex;
            mLowerHomopolymerEndIndex = lowerHpEndIndex;
            mLowerRefAdjustCount = lowerRefAdjustCount;
            mUpperHomopolymerStartIndex = lowerHpEndIndex + 1;
            mUpperHomopolymerEndIndex = upperHpEndIndex;
            mUpperRefAdjustCount = upperRefAdjustCount;
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return 0; }
    }

    private class MicrosatelliteAdjustment extends UltimaQualModel
    {
        public MicrosatelliteAdjustment()
        {
            super(UltimaVariantModel.MICROSAT_ADJUSTMENT);
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return UltimaBamUtils.ULTIMA_MAX_QUAL; }
    }

    private class SnvMnv extends UltimaQualModel
    {
        private final int mHomopolymerStartIndex;
        private final int mHomopolymerEndIndex;
        private final int mRefAdjustCount;

        public SnvMnv(
                final SimpleVariant variant, final byte[] refBases,
                final int homopolymerStartIndex, final int homopolymerEndIndex, final int refAdjustCount)
        {
            super(UltimaVariantModel.SNV);
            mHomopolymerStartIndex = homopolymerStartIndex;
            mHomopolymerEndIndex = homopolymerEndIndex;
            mRefAdjustCount = refAdjustCount;
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return 0; }

    }
}
