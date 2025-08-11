package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_TP;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_0_BOOST;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.MAX_HOMOPOLYMER;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.findHomopolymerLength;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelBuilder
{
    private final RefGenomeInterface mRefGenome;

    public static void setReadContextUltimaModels(final RefGenomeInterface refGenome, final VariantReadContext readContext)
    {
        UltimaQualModelBuilder qualModelBuilder = new UltimaQualModelBuilder(refGenome);

        byte[] coreBases = Arrays.subsetArray(
                readContext.ReadBases, readContext.VarIndex - 1, readContext.VarIndex + 1);

        UltimaQualModel qualModel = qualModelBuilder.buildContext(readContext.variant(), coreBases);

        UltimaRealignedQualModels qualModels;

        if(isCleanHpTransition(readContext, qualModel))
        {
            qualModels = new UltimaRealignedQualModels(readContext, qualModelBuilder, Lists.newArrayList());
        }
        else
        {
            qualModels = UltimaRealignedQualModelBuilder.buildUltimaRealignedQualModels(readContext, qualModelBuilder);
        }

        readContext.setUltimaRealignedQualModels(qualModels);
    }

    private static boolean isCleanHpTransition(VariantReadContext readContext, final UltimaQualModel qualModel)
    {
        if(qualModel == null || !qualModel.type().equals(UltimaModelType.HOMOPOLYMER_TRANSITION))
            return false;

        byte[] coreReadBases = Arrays.subsetArray(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        if(coreReadBases.length != readContext.RefBases.length + readContext.variant().indelLength())
            return false;

        for(int i = 0; i < coreReadBases.length; ++i)
        {
            int offset = i > readContext.leftCoreLength() ? readContext.variant().indelLength() : 0;
            if(coreReadBases[i] != readContext.RefBases[i - offset])
                return false;
        }

        return true;
    }

    public UltimaQualModelBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public UltimaQualModel buildContext(final SimpleVariant variant, final byte[] coreBases)
    {
        return buildContext(variant, coreBases, null);
    }

    public UltimaQualModel buildContext(final SimpleVariant variant, final byte[] coreBases, @Nullable final List<RefMask> refMasks)
    {
        int maxHomopolymerLength = Math.max(variant.ref().length(), MAX_HOMOPOLYMER);
        int refBaseEnd = variant.Position + maxHomopolymerLength + 1;

        // extract sufficient ref bases to set the context for most scenarios (only not for homopolymer transition)
        final byte[] origRefBases = mRefGenome.getBases(variant.Chromosome, variant.Position - 1, refBaseEnd);
        final byte[] refBases = Arrays.copyArray(origRefBases);
        if(refMasks != null)
        {
            for(RefMask refMask : refMasks)
            {
                int startIndex = Math.max(refMask.PosStart - variant.Position + 1, 0);
                int endIndex = min(refMask.PosEnd - variant.Position + 1, refBases.length - 1);
                for(int i = startIndex; i <= endIndex; i++)
                {
                    refBases[i] = refMask.BaseMask;
                }
            }
        }

        int refVarIndex = 1;

        if(variant.isIndel())
        {
            if(variant.isDelete() && isHomopolymerDeletion(variant, origRefBases))
            {
                return new UltimaHomopolymerDeletion(variant, origRefBases[2], coreBases[1], coreBases[2]);
            }
            else if(variant.isDelete() && isHomopolymerDeletion(variant, refBases))
            {
                return new UltimaHomopolymerDeletion(variant, refBases[2], coreBases[1], coreBases[2]);
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
            else
            {
                return new MicrosatelliteAdjustment();
            }
        }
        else if(variant.isSNV())
        {
            return new UltimaSnv(variant, refBases, refVarIndex, coreBases[0], coreBases[2], mRefGenome);
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

    private class MicrosatelliteAdjustment extends UltimaQualModel
    {
        public MicrosatelliteAdjustment()
        {
            super(UltimaModelType.MICROSATELLITE);
        }

        public byte calculateQual(final SAMRecord record, int varReadIndex) { return ULTIMA_MAX_QUAL_TP + ULTIMA_TP_0_BOOST; }
    }
}
