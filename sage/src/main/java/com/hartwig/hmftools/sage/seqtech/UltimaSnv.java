package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_T0;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_TP;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_0_BOOST;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.MAX_HOMOPOLYMER;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.findHomopolymerLength;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import htsjdk.samtools.SAMRecord;

class UltimaSnv extends UltimaQualModel
{
    private final UltimaHomopolymerAdjustment mLeftAdjust;
    private final UltimaHomopolymerAdjustment mRightAdjust;
    private final UltimaHomopolymerDeletion mLeftDeletion;
    private final UltimaHomopolymerDeletion mRightDeletion;

    public UltimaSnv(
            final SimpleVariant variant, final byte[] refBases, final int refVarIndex, final byte leftReadBase, final byte rightReadBase,
            final RefGenomeInterface refGenome)
    {
        super(UltimaModelType.SNV);

        byte refBase = (byte) variant.ref().charAt(0);
        byte altBase = (byte) variant.alt().charAt(0);
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
        UltimaHomopolymerAdjustment leftAdjust = null;
        UltimaHomopolymerAdjustment rightAdjust = null;
        UltimaHomopolymerDeletion leftDeletion = null;
        UltimaHomopolymerDeletion rightDeletion = null;

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
            final byte[] lowerRefBases = refGenome.getBases(variant.Chromosome, lowerRefBaseStart, lowerRefBaseEnd);
            byte hpBase = lowerRefBases[lowerRefBases.length - 1];
            lowerHpLength = findHomopolymerLength(lowerRefBases, hpBase, lowerRefBases.length - 1, false);
        }

        if(leftMatchesAlt || !leftSideDelete)
        {
            // is a HP expansion of 1 base on the left - find DP length start with the left position
            int newHpLength = lowerHpLength + 1;

            int leftHpEndIndex = 0;
            int leftHpStartIndex = leftHpEndIndex - newHpLength + 1;
            leftAdjust = new UltimaHomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, -1);
        }
        else if(leftMatchesRef)
        {
            // is a HP contraction of 1 base on the left, starting with the variant's position
            int newLowerHpLength = lowerHpLength - 1; // reflecting the new length of the HP after the delete

            int leftHpEndIndex = -1;
            int leftHpStartIndex = leftHpEndIndex - newLowerHpLength + 1;
            leftAdjust = new UltimaHomopolymerAdjustment(leftHpStartIndex, leftHpEndIndex, 1);
        }
        else
        {
            // HP deletion
            leftDeletion = new UltimaHomopolymerDeletion(variant, refBases[1], leftReadBase, rightReadBase);
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
            rightAdjust = new UltimaHomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, -1);
        }
        else if(rightMatchesRef)
        {
            // is a HP contraction of 1 base on the right, starting with the variant's position
            int newHpLength = upperHpLength - 1;

            int rightHpStartIndex = 1;
            int rightHpEndIndex = rightHpStartIndex + newHpLength - 1;
            rightAdjust = new UltimaHomopolymerAdjustment(rightHpStartIndex, rightHpEndIndex, 1);
        }
        else
        {
            rightDeletion = new UltimaHomopolymerDeletion(variant, refBases[1], leftReadBase, rightReadBase);
        }

        mLeftAdjust = leftAdjust;
        mRightAdjust = rightAdjust;
        mLeftDeletion = leftDeletion;
        mRightDeletion = rightDeletion;
    }

    public byte calculateQual(final SAMRecord record, int varReadIndex)
    {
        if(mLeftAdjust == null && mLeftDeletion == null || mRightAdjust == null && mRightDeletion == null)
        {
            return ULTIMA_MAX_QUAL_T0;
        }

        int leftQual = mLeftAdjust != null
                ?
                min(mLeftAdjust.calculateQual(record, varReadIndex), ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST)
                : min(mLeftDeletion.calculateQual(record, varReadIndex), ULTIMA_MAX_QUAL_T0);

        if(leftQual == ULTIMA_INVALID_QUAL)
        {
            return ULTIMA_INVALID_QUAL;
        }

        int rightQual = mRightAdjust != null
                ?
                min(mRightAdjust.calculateQual(record, varReadIndex), ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST)
                : min(mRightDeletion.calculateQual(record, varReadIndex), ULTIMA_MAX_QUAL_T0);

        if(rightQual == ULTIMA_INVALID_QUAL)
        {
            return ULTIMA_INVALID_QUAL;
        }

        return (byte) Math.max(leftQual, rightQual);
    }

    public String toString()
    {
        if(mLeftDeletion == null && mLeftAdjust == null)
        {
            return "non-HP";
        }

        return format("left(%s) right(%s)",
                mLeftAdjust != null ? mLeftAdjust : mLeftDeletion, mRightAdjust != null ? mRightAdjust : mRightDeletion);
    }
}
