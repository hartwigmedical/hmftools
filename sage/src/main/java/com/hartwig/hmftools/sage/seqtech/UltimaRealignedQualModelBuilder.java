package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;
import static com.hartwig.hmftools.sage.seqtech.UltimaQualModelBuilder.getStraddlingReadBases;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.INVALID_BASE;

import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class UltimaRealignedQualModelBuilder
{
    public static MergedHomopolymers mergeSandwichedHomopolymers(
            final VariantReadContext readContext, final List<Homopolymer> refHomopolymers0, final List<Homopolymer> readHomopolymers0)
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(refHomopolymers0);
        List<Homopolymer> readHomopolymers = Lists.newArrayList(readHomopolymers0);

        MergedHomopolymers mergedHomopolymers = new MergedHomopolymers();
        int refIndex = 0;
        int readIndex = 0;
        int readBasesConsumed = 0;
        int refBasesConsumed = 0;

        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);
            int refRemaining = refHomopolymers.size() - refIndex - 1;
            int readRemaining = readHomopolymers.size() - readIndex - 1;

            if(refHomopolymer.Base == readHomopolymer.Base)
            {
                int netLength = readHomopolymer.Length - refHomopolymer.Length;
                if(netLength >= 2)
                {
                    // look forward in ref
                    boolean contracted = false;
                    if(refIndex < refHomopolymers.size() - 1)
                    {
                        int extraRefIndex = 1;
                        Homopolymer forwardRefHompolymer = refHomopolymers.get(refIndex + extraRefIndex);
                        int extraRefBasesConsumed = forwardRefHompolymer.Length;

                        while(true)
                        {
                            if(forwardRefHompolymer.Base == refHomopolymer.Base)
                            {
                                int maskPositionStart = readContext.CorePositionStart + refBasesConsumed + refHomopolymer.Length;
                                int refMaskLength = IntStream.range(refIndex + 1, refIndex + extraRefIndex)
                                        .map(i -> refHomopolymers.get(i).Length)
                                        .sum();

                                int mergedLength = refMaskLength + refHomopolymers.get(refIndex).Length +
                                        refHomopolymers.get(refIndex + extraRefIndex).Length;

                                mergedHomopolymers.addRefMask(
                                        new RefMask(maskPositionStart,
                                                maskPositionStart + refMaskLength - 1, refHomopolymer.Base));

                                refHomopolymers.set(refIndex, new Homopolymer(refHomopolymer.Base, mergedLength));
                                for(int i = refIndex + 1; i <= refIndex + extraRefIndex; i++)
                                {
                                    refHomopolymers.remove(refIndex + 1);
                                }
                                contracted = true;

                                if(!readContext.variant().isIndel())
                                {
                                    int varIndexInCore = readContext.VarIndex - readContext.CoreIndexStart;
                                    if(varIndexInCore >= readBasesConsumed && varIndexInCore < readBasesConsumed + readHomopolymer.Length)
                                    {
                                        mergedHomopolymers.setVariantInMergedHomopolymers();
                                    }
                                }

                                break;
                            }
                            else if(extraRefBasesConsumed >= netLength)
                            {
                                break;
                            }

                            extraRefIndex++;
                            if(refIndex + extraRefIndex >= refHomopolymers.size())
                            {
                                break;
                            }

                            forwardRefHompolymer = refHomopolymers.get(refIndex + extraRefIndex);
                            extraRefBasesConsumed += forwardRefHompolymer.Length;
                        }
                    }

                    if(!contracted)
                    {
                        mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                        mergedHomopolymers.ReadHomopolymers.add(readHomopolymer);
                        refIndex++;
                        refBasesConsumed += refHomopolymer.Length;
                        readIndex++;
                        readBasesConsumed += readHomopolymer.Length;
                    }
                }
                else if(netLength <= -2)
                {
                    // look forward in ref
                    boolean contracted = false;
                    if(readIndex < readHomopolymers.size() - 1)
                    {
                        int extraReadIndex = 1;
                        Homopolymer forwardReadHompolymer = readHomopolymers.get(readIndex + extraReadIndex);
                        int extraReadBasesConsumed = forwardReadHompolymer.Length;

                        while(true)
                        {
                            if(forwardReadHompolymer.Base == refHomopolymer.Base)
                            {
                                int mergedLength = IntStream.range(readIndex, readIndex + extraReadIndex + 1)
                                        .map(i -> readHomopolymers.get(i).Length)
                                        .sum();

                                readHomopolymers.set(readIndex, new Homopolymer(refHomopolymer.Base, mergedLength));
                                for(int i = readIndex + 1; i <= readIndex + extraReadIndex; i++)
                                {
                                    readHomopolymers.remove(readIndex + 1);
                                }
                                contracted = true;

                                if(!readContext.variant().isIndel())
                                {
                                    int varIndexInCore = readContext.VarIndex - readContext.CoreIndexStart;
                                    if(varIndexInCore >= readBasesConsumed && varIndexInCore < readBasesConsumed + mergedLength)
                                    {
                                        mergedHomopolymers.setVariantInMergedHomopolymers();
                                    }
                                }

                                break;
                            }
                            else if(extraReadBasesConsumed >= -netLength)
                            {
                                break;
                            }

                            extraReadIndex++;
                            if(readIndex + extraReadIndex >= readHomopolymers.size())
                            {
                                break;
                            }

                            forwardReadHompolymer = readHomopolymers.get(readIndex + extraReadIndex);
                            extraReadBasesConsumed += forwardReadHompolymer.Length;
                        }
                    }

                    if(!contracted)
                    {
                        mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                        mergedHomopolymers.ReadHomopolymers.add(readHomopolymer);
                        refIndex++;
                        refBasesConsumed += refHomopolymer.Length;
                        readIndex++;
                        readBasesConsumed += readHomopolymer.Length;
                    }
                }
                else
                {
                    mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                    mergedHomopolymers.ReadHomopolymers.add(readHomopolymer);
                    refIndex++;
                    refBasesConsumed += refHomopolymer.Length;
                    readIndex++;
                    readBasesConsumed += readHomopolymer.Length;
                }
            }
            else if(readRemaining <= refRemaining)
            {
                // novel delete
                mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                refIndex++;
                refBasesConsumed += refHomopolymer.Length;
            }
            else
            {
                // novel insert
                mergedHomopolymers.ReadHomopolymers.add(readHomopolymer);
                readIndex++;
                readBasesConsumed += readHomopolymer.Length;
            }
        }

        while(refIndex < refHomopolymers.size())
        {
            mergedHomopolymers.RefHomopolymers.add(refHomopolymers.get(refIndex++));
        }

        while(readIndex < readHomopolymers.size())
        {
            mergedHomopolymers.ReadHomopolymers.add(readHomopolymers.get(readIndex++));
        }

        return mergedHomopolymers;
    }

    protected static UltimaRealignedQualModels buildUltimaRealignedQualModels(
            final VariantReadContext readContext, final UltimaQualModel baseQualModel,
            final UltimaQualModelBuilder ultimaQualModelBuilder, boolean skipSandwichMasking)
    {
        List<Homopolymer> refHomopolymers = getHomopolymers(readContext.RefBases, 0, readContext.RefBases.length - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);

        MergedHomopolymers mergedHomopolymers;
        if(skipSandwichMasking)
        {
            mergedHomopolymers = new MergedHomopolymers(refHomopolymers, readHomopolymers, false);
        }
        else
        {
            mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers);
        }

        List<UltimaRealignedQualModel> realignedQualModels = getRealignedVariants(
                readContext,
                ultimaQualModelBuilder,
                mergedHomopolymers.RefHomopolymers,
                mergedHomopolymers.ReadHomopolymers,
                mergedHomopolymers.refMasks());

        // doesn't seem a need for sorting any more, but for now leave a log warning
        realignedQualModels.sort(Comparator.comparingInt(x -> x.variant().Position));

        if(realignedQualModels.isEmpty() && !mergedHomopolymers.variantInMergedHomopolymers() && !skipSandwichMasking)
        {
            return buildUltimaRealignedQualModels(readContext, baseQualModel, ultimaQualModelBuilder, true);
        }
        else if(realignedQualModels.isEmpty() && !mergedHomopolymers.variantInMergedHomopolymers())
        {
            SG_LOGGER.info("readContext({}}) is expected to have realigned ultima variants, but none have been found", readContext);
            return new UltimaRealignedQualModels(baseQualModel);
        }

        return new UltimaRealignedQualModels(baseQualModel, realignedQualModels);
    }

    private static class HomopolymerVariant
    {
        public final byte Base;
        public final int Length;
        public final int VariantPos;
        public final int ReadCoreIndex;
        public final byte PriorBase;

        public HomopolymerVariant(byte base, int length, int variantPos, int readCoreIndex, byte priorBase)
        {
            Base = base;
            Length = length;
            VariantPos = variantPos;
            PriorBase = priorBase;
            ReadCoreIndex = readCoreIndex;
        }
    }

    private static void createRealignedQualModels(
            final List<UltimaRealignedQualModel> realignedQualModels,
            final UltimaQualModelBuilder ultimaQualModelBuilder, final VariantReadContext readContext, final List<RefMask> refMasks,
            final List<HomopolymerVariant> delHomopolymers, final List<HomopolymerVariant> insHomopolymers,
            final int lastMatchedReadCoreIndex, final int lastMatchedRefPos, final byte lastMatchedBase, final int lastMatchedLength,
            final byte priorRefBase, final byte priorReadBase)
    {
        if(delHomopolymers.isEmpty() && insHomopolymers.isEmpty())
            return;

        if(lastMatchedRefPos == -1 || lastMatchedReadCoreIndex == -1 || lastMatchedLength == -1)
        {
            throw new IllegalArgumentException(
                    "first core homopolymers of the read and ref do not have the same base, variant: " + readContext.variant().toString());
        }

        String chromosome = readContext.variant().chromosome();
        int origReadCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;
        for(int i = 0; i < delHomopolymers.size(); i++)
        {
            HomopolymerVariant delHomopolymer = delHomopolymers.get(i);

            int varReadIndexOffset = lastMatchedReadCoreIndex - origReadCoreIndex;
            String delBasesString = String.valueOf((char) delHomopolymer.Base).repeat(delHomopolymer.Length);
            byte firstRefBase = delHomopolymer.PriorBase;
            // in GGGGAG -> GGGG, don't want to attach ins G to start of poly-G
            byte firstAltBase = (delHomopolymer.Base == lastMatchedBase && i == 0) ? priorReadBase : lastMatchedBase;
            if(firstRefBase == INVALID_BASE || firstAltBase == INVALID_BASE)
            {
                firstRefBase = readContext.readBasesBytes()[readContext.CoreIndexStart - 1];
                firstAltBase = firstRefBase;
            }

            if(delHomopolymer.Base == lastMatchedBase && i == 0)  // in GGGGAG -> GGGG, don't want to attach del G to start of poly-G
            {
                varReadIndexOffset -= lastMatchedLength;
            }

            String ref = String.valueOf((char)firstRefBase) + delBasesString;
            String alt = String.valueOf((char)firstAltBase);

            SimpleVariant variant = new SimpleVariant(chromosome, delHomopolymer.VariantPos, ref, alt);

            int adjustedVarReadIndex = readContext.VarIndex + varReadIndexOffset;
            byte[] straddlingBases = getStraddlingReadBases(variant, readContext.readBasesBytes(), adjustedVarReadIndex);

            int varIndex = varReadIndexOffset + readContext.VarIndex;
            UltimaQualModel baseQualModel = ultimaQualModelBuilder.buildContext(variant, straddlingBases, refMasks);

            UltimaRealignedQualModel realignedQualModel = new UltimaRealignedQualModel(
                    variant, baseQualModel, varReadIndexOffset, varIndex,
                    delHomopolymer.VariantPos - readContext.CorePositionStart);

            realignedQualModels.add(realignedQualModel);
        }

        for(int i = 0; i < insHomopolymers.size(); i++)
        {
            HomopolymerVariant insHomopolymer = insHomopolymers.get(i);

            int varReadIndexOffset = insHomopolymer.ReadCoreIndex - origReadCoreIndex;
            String insBasesString = String.valueOf((char) insHomopolymer.Base).repeat(insHomopolymer.Length);
            // in GGGG -> GGGGAG, don't want to attach ins G to start of poly-G
            byte firstRefBase = (insHomopolymer.Base == lastMatchedBase && i == 0) ? priorRefBase : lastMatchedBase;
            byte firstAltBase = insHomopolymer.PriorBase;
            if(firstRefBase == INVALID_BASE || firstAltBase == INVALID_BASE)
            {
                firstRefBase = readContext.readBasesBytes()[readContext.CoreIndexStart - 1];
                firstAltBase = firstRefBase;
            }

            String ref = String.valueOf((char)firstRefBase);
            String alt = (char)firstAltBase + insBasesString;

            SimpleVariant variant = new SimpleVariant(chromosome, insHomopolymer.VariantPos, ref, alt);

            int adjustedVarReadIndex = readContext.VarIndex + varReadIndexOffset;
            byte[] straddlingBases = getStraddlingReadBases(variant, readContext.readBasesBytes(), varReadIndexOffset);

            UltimaQualModel baseQualModel = ultimaQualModelBuilder.buildContext(variant, straddlingBases, refMasks);

            UltimaRealignedQualModel realignedQualModel = new UltimaRealignedQualModel(
                    variant, baseQualModel, varReadIndexOffset, adjustedVarReadIndex,
                    insHomopolymer.VariantPos - readContext.CorePositionStart);

            realignedQualModels.add(realignedQualModel);
        }
    }

    private static void extendHomopolymerVariants(final List<HomopolymerVariant> variants, HomopolymerVariant variant)
    {
        if(variants.isEmpty())
        {
            variants.add(variant);
            return;
        }

        HomopolymerVariant lastVariant = variants.get(variants.size() - 1);
        if(lastVariant.Base == variant.Base)
        {
            variants.set(
                    variants.size() - 1,
                    new HomopolymerVariant(
                            lastVariant.Base,
                            lastVariant.Length + variant.Length,
                            lastVariant.VariantPos,
                            lastVariant.ReadCoreIndex,
                            lastVariant.PriorBase));

            return;
        }

        variants.add(variant);
    }

    @VisibleForTesting
    public static List<UltimaRealignedQualModel> getRealignedVariants(final VariantReadContext readContext,
            final UltimaQualModelBuilder ultimaQualModelBuilder, final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers, final List<RefMask> refMasks)
    {
        List<UltimaRealignedQualModel> realignedVariants = Lists.newArrayList();

        int lastMatchedReadCoreIndex = -1;
        int readCoreIndex = 0;
        byte lastMatchedBase = 0;
        int lastMatchedRefPos = -1;
        int lastMatchedLength = -1;
        byte priorRefBase = -1;
        byte priorReadBase = -1;
        int refIndex = 0;
        int readIndex = 0;
        int refPos = readContext.CorePositionStart;
        List<HomopolymerVariant> delHomopolymers = Lists.newArrayList();
        List<HomopolymerVariant> insHomopolymers = Lists.newArrayList();
        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolymer.Base && refHomopolymer.Length == readHomopolymer.Length)
            {
                createRealignedQualModels(realignedVariants, ultimaQualModelBuilder, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

                lastMatchedLength = refHomopolymer.Length;
                priorRefBase = refIndex == 0 ? -1 : refHomopolymers.get(refIndex - 1).Base;
                priorReadBase = readIndex == 0 ? -1 : readHomopolymers.get(readIndex - 1).Base;

                lastMatchedBase = refHomopolymer.Base;
                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                lastMatchedRefPos = refPos - 1;
                readCoreIndex += readHomopolymer.Length;
                lastMatchedReadCoreIndex = readCoreIndex - 1;

                delHomopolymers.clear();
                insHomopolymers.clear();
                continue;
            }

            if(refHomopolymer.Base == readHomopolymer.Base)
            {
                createRealignedQualModels(realignedVariants, ultimaQualModelBuilder, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

                lastMatchedLength = min(refHomopolymer.Length, readHomopolymer.Length);
                priorRefBase = refIndex == 0 ? -1 : refHomopolymers.get(refIndex - 1).Base;
                priorReadBase = readIndex == 0 ? -1 : readHomopolymers.get(readIndex - 1).Base;

                lastMatchedBase = refHomopolymer.Base;
                ++refIndex;
                ++readIndex;
                lastMatchedRefPos = refPos + min(refHomopolymer.Length, readHomopolymer.Length) - 1;
                refPos += refHomopolymer.Length;
                lastMatchedReadCoreIndex = readCoreIndex + min(refHomopolymer.Length, readHomopolymer.Length) - 1;
                readCoreIndex += readHomopolymer.Length;

                delHomopolymers.clear();
                insHomopolymers.clear();

                if(refHomopolymer.Length < readHomopolymer.Length)
                {
                    byte priorBase = readIndex <= 1 ? INVALID_BASE : readHomopolymers.get(readIndex - 2).Base;
                    HomopolymerVariant insHomopolymer = new HomopolymerVariant(
                            readHomopolymer.Base,
                            readHomopolymer.Length - refHomopolymer.Length,
                            refPos - refHomopolymer.Length - 1,
                            readCoreIndex - readHomopolymer.Length - 1,
                            priorBase);

                    insHomopolymers.add(insHomopolymer);
                }
                else
                {
                    byte priorBase = refIndex <= 1 ? INVALID_BASE : refHomopolymers.get(refIndex - 2).Base;
                    HomopolymerVariant delHomopolymer = new HomopolymerVariant(
                            refHomopolymer.Base,
                            refHomopolymer.Length - readHomopolymer.Length,
                            refPos - refHomopolymer.Length - 1,
                            readCoreIndex - readHomopolymer.Length - 1,
                            priorBase);

                    delHomopolymers.add(delHomopolymer);
                }

                continue;
            }

            int refHomopolymersLeft = refHomopolymers.size() - refIndex - 1;
            int readHomopolymersLeft = readHomopolymers.size() - readIndex - 1;
            if(refHomopolymersLeft == readHomopolymersLeft)
            {
                byte priorBase = refIndex == 0 ? INVALID_BASE : refHomopolymers.get(refIndex - 1).Base;
                HomopolymerVariant delHomopolymer = new HomopolymerVariant(
                        refHomopolymer.Base, refHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

                extendHomopolymerVariants(delHomopolymers, delHomopolymer);

                priorBase = readIndex == 0 ? INVALID_BASE : readHomopolymers.get(readIndex - 1).Base;
                HomopolymerVariant insHomopolyer = new HomopolymerVariant(
                        readHomopolymer.Base, readHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

                extendHomopolymerVariants(insHomopolymers, insHomopolyer);

                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                readCoreIndex += readHomopolymer.Length;
                continue;
            }

            if(refHomopolymersLeft > readHomopolymersLeft)
            {
                byte priorBase = refIndex == 0 ? INVALID_BASE : refHomopolymers.get(refIndex - 1).Base;
                HomopolymerVariant delHomopolymer = new HomopolymerVariant(
                        refHomopolymer.Base, refHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

                extendHomopolymerVariants(delHomopolymers, delHomopolymer);
                ++refIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            byte priorBase = readIndex == 0 ? INVALID_BASE : readHomopolymers.get(readIndex - 1).Base;
            HomopolymerVariant insHomopolymer = new HomopolymerVariant(
                    readHomopolymer.Base, readHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

            extendHomopolymerVariants(insHomopolymers, insHomopolymer);
            ++readIndex;
            readCoreIndex += readHomopolymer.Length;
        }

        while(readIndex < readHomopolymers.size())
        {
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);
            byte priorBase = readIndex == 0 ? INVALID_BASE : readHomopolymers.get(readIndex - 1).Base;
            HomopolymerVariant insHomopolymer = new HomopolymerVariant(
                    readHomopolymer.Base, readHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

            extendHomopolymerVariants(insHomopolymers, insHomopolymer);
            ++readIndex;
            readCoreIndex += readHomopolymer.Length;
        }

        while(refIndex < refHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            byte priorBase = refIndex == 0 ? INVALID_BASE : refHomopolymers.get(refIndex - 1).Base;
            HomopolymerVariant delHomopolymer = new HomopolymerVariant(
                    refHomopolymer.Base, refHomopolymer.Length, refPos - 1, readCoreIndex - 1, priorBase);

            extendHomopolymerVariants(delHomopolymers, delHomopolymer);
            ++refIndex;
            refPos += refHomopolymer.Length;
        }

        createRealignedQualModels(realignedVariants, ultimaQualModelBuilder, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

        return realignedVariants;
    }
}
