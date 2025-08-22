package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.sage.seqtech.UltimaRealignedQualModelBuilder.buildUltimaRealignedQualModels;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.MAX_HOMOPOLYMER;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.findHomopolymerLength;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelBuilder
{
    private final RefGenomeInterface mRefGenome;

    public static void setReadContextUltimaModels(final RefGenomeInterface refGenome, final ReadContextCounter readContextCounter)
    {
        VariantReadContext readContext = readContextCounter.readContext();
        UltimaQualModelBuilder qualModelBuilder = new UltimaQualModelBuilder(refGenome);

        byte[] straddlingReadBases = Arrays.subsetArray(
                readContext.ReadBases, readContext.VarIndex - 1, readContext.VarIndex + 1);

        UltimaQualModel qualModel = qualModelBuilder.buildContext(readContext.variant(), straddlingReadBases);

        UltimaRealignedQualModels qualModels;

        if(!canSkipRealignedModels(readContext))
        {
            qualModels = new UltimaRealignedQualModels(qualModel);
        }
        else
        {
            qualModels = buildUltimaRealignedQualModels(readContext, qualModel, qualModelBuilder, false);;
        }

        readContextCounter.ultimaData().setQualModels(qualModels);
    }

    private static boolean canSkipRealignedModels(final VariantReadContext readContext)
    {
        int indelLength = readContext.variant().indelLength();
        int coreLength = readContext.coreLength();

        if(coreLength != readContext.RefBases.length + indelLength)
            return false;

        int readIndex = readContext.CoreIndexStart;

        if(readContext.variant().isIndel())
        {
            int refIndex = 0;

            while(readIndex < readContext.ReadBases.length && refIndex < readContext.RefBases.length)
            {
                if(readIndex == readContext.VarIndex + 1)
                {
                    if(readContext.variant().isDelete())
                        refIndex += abs(indelLength);
                    else
                        readIndex += indelLength;
                }

                if(readContext.ReadBases[readIndex] != readContext.RefBases[refIndex])
                    return false;

                ++readIndex;
                ++refIndex;
            }
        }
        else
        {
            int altLength = readContext.variant().altLength();

            for(int i = 0; i < coreLength; ++i, ++readIndex)
            {
                if(readIndex >= readContext.VarIndex && readIndex <= readContext.VarIndex + (altLength - 1))
                    continue;

                if(readContext.ReadBases[readIndex] != readContext.RefBases[i])
                    return false;
            }
        }

        return true;
    }

    public UltimaQualModelBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public UltimaQualModel buildContext(final SimpleVariant variant, final byte[] straddlingReadBases)
    {
        return buildContext(variant, straddlingReadBases, null);
    }

    public UltimaQualModel buildContext(final SimpleVariant variant, final byte[] straddlingReadBases, @Nullable final List<RefMask> refMasks)
    {
        // straddling read bases are the 3 read bases around the variant itself
        int maxHomopolymerLength = max(variant.ref().length(), MAX_HOMOPOLYMER);
        int refBaseEnd = variant.Position + maxHomopolymerLength + 1;

        // extract sufficient ref bases to set the context for most scenarios (only not for homopolymer transition)

        // TODO: could use the variant read context's ref bases if long enough, rather than look up the ref genome again

        byte[] origRefBases = mRefGenome.getBases(variant.Chromosome, variant.Position - 1, refBaseEnd);

        byte[] refBases;
        if(refMasks != null)
        {
            refBases = Arrays.copyArray(origRefBases);

            for(RefMask refMask : refMasks)
            {
                int startIndex = max(refMask.PosStart - variant.Position + 1, 0);
                int endIndex = min(refMask.PosEnd - variant.Position + 1, refBases.length - 1);
                for(int i = startIndex; i <= endIndex; i++)
                {
                    refBases[i] = refMask.BaseMask;
                }
            }
        }
        else
        {
            refBases = origRefBases;
        }

        int refVarIndex = 1; // since the ref bases above start 1 base before the variant

        if(variant.isIndel())
        {
            if(variant.isDelete() && isHomopolymerDeletion(variant, origRefBases))
            {
                // HP base is the first ref base after the variant's position, ie the first base deleted
                // the straddling read bases here are the variant's ref and first ref after the delete
                return UltimaHomopolymerDeletion.fromDelete(origRefBases[2], straddlingReadBases[1], straddlingReadBases[2]);
                // return new UltimaHomopolymerDeletion(variant, origRefBases[2], straddlingReadBases[1], straddlingReadBases[2]);
            }
            else if(variant.isDelete() && isHomopolymerDeletion(variant, refBases))
            {
                return UltimaHomopolymerDeletion.fromDelete(refBases[2], straddlingReadBases[1], straddlingReadBases[2]);
                // return new UltimaHomopolymerDeletion(variant, refBases[2], straddlingReadBases[1], straddlingReadBases[2]);
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
                    return new UltimaHomopolymerTransitionDeletion(variant, homopolymerTransitionIndex, lowerHpLength, upperHpLength);
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

                return new UltimaHomopolymerAdjustment(homopolymerStartIndex, homopolymerEndIndex, refAdjustCount);
            }
        }
        else if(variant.isSNV())
        {
            UltimaSnv snvModel = new UltimaSnv(variant, refBases, refVarIndex, straddlingReadBases[0], straddlingReadBases[2], mRefGenome);

            if(snvModel.canCompute())
                return snvModel;
        }
        else
        {
            // MNV case which may be a base shift
            UltimaBaseShift mnvModel = new UltimaBaseShift(variant, refBases, refVarIndex, straddlingReadBases[0], straddlingReadBases[2], mRefGenome);

            if(mnvModel.canCompute())
                return mnvModel;
        }

        // return a fall-back for no valid model
        return new OtherVariant();
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

    private class OtherVariant extends UltimaQualModel
    {
        public OtherVariant()
        {
            super(UltimaModelType.OTHER);
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return INVALID_BASE_QUAL; }

        @Override
        public boolean canCompute() { return false; }
    }
}
