package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.BQR_CACHE;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.MAX_HOMOPOLYMER;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.findHomopolymerLength;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import htsjdk.samtools.SAMRecord;

class UltimaBaseShift extends UltimaQualModel
{
    private final int mMnvLength;
    private final UltimaHomopolymerAdjustment mLeftAdjust;
    private final UltimaHomopolymerAdjustment mRightAdjust;
    private final UltimaHomopolymerDeletion mLeftDeletion;
    private final UltimaHomopolymerDeletion mRightDeletion;

    public UltimaBaseShift(
            final SimpleVariant variant, final byte[] refBases, final int refVarIndex, final byte leftReadBase, final byte rightReadBase,
            final RefGenomeInterface refGenome)
    {
        super(UltimaModelType.BASE_SHIFT);

        // check bases on either side of the variant
        int mnvLength = variant.refLength();
        mMnvLength = mnvLength;

        // conditions for a base shift is a deletion or contraction on one side, remaining bases shifting and last MNV base an insert

        // test each side in turn
        boolean validShiftLeft = true;

        for(int i = 0; i < mnvLength - 1; ++i)
        {
            // ie ACG -> CGT, or AC -> CT
            char refBase = variant.ref().charAt(i + 1);
            char altBase = variant.alt().charAt(i);

            if(refBase != altBase)
            {
                validShiftLeft = false;
                break;
            }
        }

        boolean validShiftRight = true;

        if(validShiftLeft)
        {
            validShiftRight = false;
        }
        else
        {
            for(int i = 0; i < mnvLength - 1; ++i)
            {
                // ie ACG -> TAC, or AC -> TA
                char refBase = variant.ref().charAt(i);
                char altBase = variant.alt().charAt(i + 1);

                if(refBase != altBase)
                {
                    validShiftRight = false;
                    break;
                }
            }
        }

        if(!validShiftLeft && !validShiftRight)
        {
            mLeftAdjust = null;
            mRightAdjust = null;
            mLeftDeletion = null;
            mRightDeletion = null;
            return;
        }

        // one side will be a delete (either partial or complete), the other a 1-base insert
        UltimaHomopolymerAdjustment leftAdjust = null;
        UltimaHomopolymerAdjustment rightAdjust = null;
        UltimaHomopolymerDeletion leftDeletion = null;
        UltimaHomopolymerDeletion rightDeletion = null;

        byte leftRefBase = refBases[refVarIndex];
        byte rightRefBase = refBases[refVarIndex + mnvLength - 1];
        byte leftAltBase = (byte)variant.alt().charAt(0);
        byte rightAltBase = (byte)variant.alt().charAt(mnvLength - 1);

        byte nextLeftRefBase = refBases[refVarIndex - 1];
        byte nextRightRefBase = refBases[refVarIndex + mnvLength];

        if(validShiftLeft)
        {
            // find the deletion or contraction on the left
            if(leftRefBase == nextLeftRefBase) // a contraction
            {
                int lowerRefBaseEnd = variant.Position;
                int lowerRefBaseStart = lowerRefBaseEnd - MAX_HOMOPOLYMER;
                final byte[] lowerRefBases = refGenome.getBases(variant.Chromosome, lowerRefBaseStart, lowerRefBaseEnd);
                int lowerHpLength = findHomopolymerLength(lowerRefBases, leftRefBase, lowerRefBases.length - 1, false);
                int newLowerHpLength = lowerHpLength - 1;

                int leftHpEndIndex = -1;
                int leftHpStartIndex = leftHpEndIndex - newLowerHpLength + 1;
                leftAdjust = new UltimaHomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, 1);
            }
            else
            {
                byte rightBase = (byte)variant.alt().charAt(1);
                leftDeletion = UltimaHomopolymerDeletion.fromMnv(leftRefBase, leftAltBase, leftReadBase, rightBase);
            }

            // check if the insertion on the right is an expansion or a novel insertion
            int rightHpStartIndex = mnvLength - 1;

            if(rightAltBase == nextRightRefBase)
            {
                // the alt has expanded an existing homopolymer
                int upperRefBaseStart = variant.positionEnd() + 1;
                int upperRefBaseEnd = upperRefBaseStart + MAX_HOMOPOLYMER;
                final byte[] upperRefBases = refGenome.getBases(variant.Chromosome, upperRefBaseStart, upperRefBaseEnd);
                int upperHpLength = findHomopolymerLength(upperRefBases, nextRightRefBase, 0, true);
                int newUpperHpLength = upperHpLength + 1;

                int rightHpEndIndex = rightHpStartIndex + newUpperHpLength - 1;
                rightAdjust = new UltimaHomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, -1);
            }
            else
            {
                // a novel insertion so is a HP expansion owf 1 base on the right
                int rightHpEndIndex = rightHpStartIndex;
                rightAdjust = new UltimaHomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, -1);
            }
        }
        else
        {
            if(rightRefBase == nextRightRefBase)
            {
                int upperRefBaseStart = variant.positionEnd();
                int upperRefBaseEnd = upperRefBaseStart + MAX_HOMOPOLYMER;
                final byte[] uppperRefBases = refGenome.getBases(variant.Chromosome, upperRefBaseStart, upperRefBaseEnd);
                int upperHpLength = findHomopolymerLength(uppperRefBases, rightRefBase, 0, true);
                int newUpperHpLength = upperHpLength - 1;

                int rightHpStartIndex = mMnvLength;
                int rightHpEndIndex = rightHpStartIndex + newUpperHpLength - 1;
                rightAdjust = new UltimaHomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, 1);
            }
            else
            {
                byte leftBase = (byte)variant.alt().charAt(mnvLength - 2);
                rightDeletion = UltimaHomopolymerDeletion.fromMnv(rightRefBase, rightAltBase, leftBase, rightReadBase);
            }

            int leftHpEndIndex = 0;

            if(leftAltBase == nextLeftRefBase)
            {
                int lowerRefBaseEnd = variant.Position - 1;
                int upperRefBaseStart = lowerRefBaseEnd - MAX_HOMOPOLYMER;
                final byte[] lowerRefBases = refGenome.getBases(variant.Chromosome, upperRefBaseStart, lowerRefBaseEnd);
                int lowerHpLength = findHomopolymerLength(lowerRefBases, nextLeftRefBase, lowerRefBases.length - 1, false);
                int newUpperHpLength = lowerHpLength + 1;

                int leftHpStartIndex = leftHpEndIndex - newUpperHpLength + 1;
                leftAdjust = new UltimaHomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, -1);
            }
            else
            {
                int leftHpStartIndex = leftHpEndIndex;
                leftAdjust = new UltimaHomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, -1);
            }
        }

        mLeftAdjust = leftAdjust;
        mRightAdjust = rightAdjust;
        mLeftDeletion = leftDeletion;
        mRightDeletion = rightDeletion;
    }

    @Override
    public boolean canCompute()
    {
        return (mLeftAdjust != null || mLeftDeletion != null) && (mRightAdjust != null || mRightDeletion != null);
    }

    public byte calculateQual(final SAMRecord record, int varReadIndex)
    {
        if(!canCompute())
            return ULTIMA_INVALID_QUAL;

        int leftQual = 0;

        if(mLeftAdjust != null)
        {
            leftQual = mLeftAdjust.calculateQual(record, varReadIndex);
        }
        else
        {
            leftQual = mLeftDeletion.calculateQual(record, varReadIndex);
        }

        if(leftQual == ULTIMA_INVALID_QUAL)
            return ULTIMA_INVALID_QUAL;

        int rightQual = 0;
        int mnvEndReadIndex = varReadIndex + mMnvLength - 1;

        if(mRightAdjust != null)
        {
            rightQual = mRightAdjust.calculateQual(record, mnvEndReadIndex);
        }
        else
        {
            rightQual = mRightDeletion.calculateQual(record, mnvEndReadIndex);
        }

        if(rightQual == ULTIMA_INVALID_QUAL)
            return ULTIMA_INVALID_QUAL;

        return (byte)max(leftQual, rightQual);
    }

    @VisibleForTesting
    public UltimaHomopolymerAdjustment leftAdjust() { return mLeftAdjust; }
    public UltimaHomopolymerAdjustment rightAdjust() { return mRightAdjust; }
    public UltimaHomopolymerDeletion leftDeletion() { return mLeftDeletion; }
    public UltimaHomopolymerDeletion rightDeletion() { return mRightDeletion; }

    public boolean isLeftShift() { return mLeftDeletion != null || mLeftAdjust.refAdjustCount() > 0; }
    public boolean isRightShift() { return mRightDeletion != null || mRightAdjust.refAdjustCount() > 0; }

    public String toString()
    {
        if(!canCompute())
            return "invalid";

        return format("mnvLength(%d) left(%s) right(%s)",
                mMnvLength, mLeftAdjust != null ? mLeftAdjust : mLeftDeletion, mRightAdjust != null ? mRightAdjust : mRightDeletion);
    }
}
