// TODO: REVIEW
package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.CYCLE_BASE_INDEX;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.INVALID_BASE;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.isCleanSnv;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.getHomopolymers;

import static htsjdk.samtools.CigarOperator.I;

import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.Homopolymer;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.MergedHomopolymers;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.RefMask;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;

public class UltimaRealignedQualModelsBuilder_0
{
    private static int cycleCount(final List<Homopolymer> homopolymers, int startIndex)
    {
        if(startIndex < 0 || startIndex >= homopolymers.size())
        {
            return 0;
        }

        Homopolymer prev_homopolymer = homopolymers.get(startIndex);
        int count = 1;
        for(int i = startIndex + 1; i < homopolymers.size(); i++)
        {
            Homopolymer homopolymer = homopolymers.get(i);
            if(CYCLE_BASE_INDEX.get(homopolymer.Base) < CYCLE_BASE_INDEX.get(prev_homopolymer.Base))
            {
                count++;
            }

            prev_homopolymer = homopolymer;
        }

        return count;
    }

    // TODO: Do we need the requireCycleShift option any more?
    @VisibleForTesting
    public static MergedHomopolymers mergeSandwichedHomopolymers(@Nullable final VariantReadContext readContext,
            final List<Homopolymer> refHomopolymers0, final List<Homopolymer> readHomopolymers0, boolean requireCycleShift)
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
                        // TODO: ASK THOMAS. Ref cycle count?
                        int extraRefIndex = 1;
                        Homopolymer forwardRefHompolymer = refHomopolymers.get(refIndex + extraRefIndex);
                        int extraRefBasesConsumed = forwardRefHompolymer.Length;
                        int refCycleCount = cycleCount(refHomopolymers, refIndex);
                        int readCycleCount = cycleCount(readHomopolymers, readIndex);
                        while(true)
                        {
                            if(forwardRefHompolymer.Base == refHomopolymer.Base && (!requireCycleShift || refCycleCount != readCycleCount))
                            {
                                int maskPositionStart = readContext.corePositionStart() + refBasesConsumed + refHomopolymer.Length;
                                int refMaskLength = IntStream.range(refIndex + 1, refIndex + extraRefIndex)
                                        .map(i -> refHomopolymers.get(i).Length)
                                        .sum();

                                int mergedLength = refMaskLength + refHomopolymers.get(refIndex).Length +
                                        refHomopolymers.get(refIndex + extraRefIndex).Length;

                                mergedHomopolymers.addRefMask(
                                        new UltimaRealignedQualModelsBuilder.RefMask(maskPositionStart,
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
                        int refCycleCount = cycleCount(refHomopolymers, refIndex);
                        int readCycleCount = cycleCount(readHomopolymers, readIndex);
                        while(true)
                        {
                            if(forwardReadHompolymer.Base == refHomopolymer.Base && (!requireCycleShift || refCycleCount != readCycleCount))
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

    public static UltimaRealignedQualModels buildUltimaRealignedQualModels(final VariantReadContext readContext,
            final UltimaQualCalculator ultimaQualCalculator)
    {
        return buildUltimaRealignedQualModels(readContext, ultimaQualCalculator, false);
    }

    private static UltimaRealignedQualModels buildUltimaRealignedQualModels(final VariantReadContext readContext,
            final UltimaQualCalculator ultimaQualCalculator, boolean skipSandwichMasking)
    {
        if(!skipSandwichMasking && isCleanSnv(readContext))
        {
            return new UltimaRealignedQualModels(readContext, ultimaQualCalculator);  // if a clean SNV, want to take max of quals, not min
        }

        if(!skipSandwichMasking && readContext.variant().isInsert())
        {
            CigarElement insertEl = new CigarElement(readContext.alt().length() - 1, I);
            List<CigarElement> coreCigarElements = readContext.coreCigarElements();
            if(!coreCigarElements.contains(insertEl))
            {
                return new UltimaRealignedQualModels(readContext, ultimaQualCalculator);
            }
        }

        List<Homopolymer> refHomopolymers = getHomopolymers(readContext.RefBases, 0, readContext.RefBases.length - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        MergedHomopolymers mergedHomopolymers;
        if(skipSandwichMasking)
        {
            mergedHomopolymers = new MergedHomopolymers(refHomopolymers, readHomopolymers, false);
        }
        else
        {
            mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers, false);
        }

        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext,
                ultimaQualCalculator,
                mergedHomopolymers.RefHomopolymers,
                mergedHomopolymers.ReadHomopolymers,
                mergedHomopolymers.refMasks());

        // TODO: Look into why we need to do this.
        realignedVariants.sort(Comparator.comparingInt(x -> x.variant().Position));

        // TODO: Move this into a unit test.
        for(int i = 1; i < realignedVariants.size(); i++)
        {
            if(realignedVariants.get(i).variant().Position < realignedVariants.get(i - 1).variant().Position)
            {
                throw new IllegalStateException(format("Realigned variants are out of order: readContext(%s)", readContext));
            }
        }

        // TODO: This has been duct taped.
        List<UltimaRealignedQualModel> realignedQualModels = realignedVariants;

        // TODO: Indel balance?
        // getQualVariants(mergedHomopolymers.variantInMergedHomopolymers(), readContext.variant(), realignedVariants);

        if(realignedQualModels.isEmpty() && !mergedHomopolymers.variantInMergedHomopolymers() && !skipSandwichMasking)
        {
            return buildUltimaRealignedQualModels(readContext, ultimaQualCalculator, true);  // TODO: should handle skipping first/last sandwiched?
        }
        else if(realignedQualModels.isEmpty() && !mergedHomopolymers.variantInMergedHomopolymers())
        {
            SG_LOGGER.info(format("readContext(%s) is expected to have realigned ultima variants, but none have been found", readContext.toString()));
            return new UltimaRealignedQualModels(readContext, ultimaQualCalculator);
        }

        return new UltimaRealignedQualModels(readContext, ultimaQualCalculator, realignedQualModels);
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

    private static void createRealignedVariants(final List<UltimaRealignedQualModel> realignedVariants,
            @Nullable final UltimaQualCalculator ultimaQualCalculator, final VariantReadContext readContext, final List<RefMask> refMasks,
            final List<HomopolymerVariant> delHomopolymers, final List<HomopolymerVariant> insHomopolymers,
            final int lastMatchedReadCoreIndex, final int lastMatchedRefPos, final byte lastMatchedBase, final int lastMatchedLength,
            final byte priorRefBase, final byte priorReadBase)
    {
        // we allow a null ultimaQualCalculator for testing only purposes

        if(delHomopolymers.isEmpty() && insHomopolymers.isEmpty())
        {
            return;
        }

        if(lastMatchedRefPos == -1 || lastMatchedReadCoreIndex == -1 || lastMatchedLength == -1)
        {
            throw new IllegalArgumentException("Looks as though the first core homopolymers of the read and ref do not have the same base, variant: " + readContext.variant().toString());
        }

        String chromosome = readContext.variant().chromosome();
        int origReadCoreIndex = readContext.varIndex() - readContext.coreIndexStart();
        for(int i = 0; i < delHomopolymers.size(); i++)
        {
            HomopolymerVariant delHomopolymer = delHomopolymers.get(i);

            int varReadIndexOffset = lastMatchedReadCoreIndex - origReadCoreIndex;
            String delBasesString = String.valueOf((char) delHomopolymer.Base).repeat(delHomopolymer.Length);
            byte firstRefBase = delHomopolymer.PriorBase;
            byte firstAltBase = delHomopolymer.Base == lastMatchedBase ? priorReadBase : lastMatchedBase;
            if(firstRefBase == INVALID_BASE || firstAltBase == INVALID_BASE)
            {
                firstRefBase = readContext.readBasesBytes()[readContext.coreIndexStart() - 1];
                firstAltBase = firstRefBase;
            }

            if(delHomopolymer.Base == lastMatchedBase && i == 0)  // in GGGGAG -> GGGG, don't want to attach del G to start of poly-G
            {
                varReadIndexOffset -= lastMatchedLength;
            }

            String ref = String.valueOf((char) firstRefBase) + delBasesString;
            String alt = String.valueOf((char) firstAltBase);

            SimpleVariant variant = new SimpleVariant(chromosome, delHomopolymer.VariantPos, ref, alt);
            byte[] coreBases = Arrays.subsetArray(
                    readContext.readBasesBytes(),
                    readContext.varIndex() + varReadIndexOffset - 1,
                    readContext.varIndex() + varReadIndexOffset + 1);

            // common scenario, we want to pass in the base after the variant, not the SNV base itself, for the right straddle base
            if(varReadIndexOffset == -1 && !readContext.variant().isIndel())
            {
                coreBases[2] = readContext.readBasesBytes()[readContext.varIndex() + varReadIndexOffset + 2];
            }

            int varIndex = varReadIndexOffset + readContext.varIndex();
            UltimaQualModel baseQualModel = ultimaQualCalculator == null
                    ? null
                    : ultimaQualCalculator.buildContext(variant, coreBases, refMasks);

            UltimaRealignedQualModel realignedQualModel = baseQualModel == null
                    ? new UltimaRealignedQualModel(variant, varReadIndexOffset)
                    : new UltimaRealignedQualModel(variant, baseQualModel, varReadIndexOffset, varIndex,
                            delHomopolymer.VariantPos - readContext.corePositionStart());

            realignedVariants.add(realignedQualModel);
        }

        for(int i = 0; i < insHomopolymers.size(); i++)
        {
            HomopolymerVariant insHomopolymer = insHomopolymers.get(i);

            int varReadIndexOffset = insHomopolymer.ReadCoreIndex - origReadCoreIndex + i;
            String insBasesString = String.valueOf((char) insHomopolymer.Base).repeat(insHomopolymer.Length);
            byte firstRefBase = insHomopolymer.Base == lastMatchedBase ? priorRefBase : lastMatchedBase;
            byte firstAltBase = insHomopolymer.PriorBase;
            if(firstRefBase == INVALID_BASE || firstAltBase == INVALID_BASE)
            {
                firstRefBase = readContext.readBasesBytes()[readContext.coreIndexStart() - 1];
                firstAltBase = firstRefBase;
            }

            if(insHomopolymer.Base == lastMatchedBase && i == 0)  // in GGGG -> GGGGAG, don't want to attach ins G to start of poly-G
            {
                varReadIndexOffset -= lastMatchedLength;
            }

            String ref = String.valueOf((char) firstRefBase);
            String alt = String.valueOf((char) firstAltBase) + insBasesString;

            SimpleVariant variant = new SimpleVariant(chromosome, insHomopolymer.VariantPos, ref, alt);
            byte[] coreBases = Arrays.subsetArray(
                    readContext.readBasesBytes(),
                    readContext.varIndex() + varReadIndexOffset - 1,
                    readContext.varIndex() + varReadIndexOffset + 1);

            int varIndex = varReadIndexOffset + readContext.varIndex();
            UltimaQualModel baseQualModel = ultimaQualCalculator == null
                    ? null
                    : ultimaQualCalculator.buildContext(variant, coreBases, refMasks);

            UltimaRealignedQualModel realignedQualModel = baseQualModel == null
                    ? new UltimaRealignedQualModel(variant, varReadIndexOffset)
                    : new UltimaRealignedQualModel(variant, baseQualModel, varReadIndexOffset, varIndex,
                            insHomopolymer.VariantPos - readContext.corePositionStart());

            realignedVariants.add(realignedQualModel);
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
            final UltimaQualCalculator ultimaQualCalculator, final List<Homopolymer> refHomopolymers,
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
        int refPos = readContext.corePositionStart();
        List<HomopolymerVariant> delHomopolymers = Lists.newArrayList();
        List<HomopolymerVariant> insHomopolymers = Lists.newArrayList();
        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolymer.Base && refHomopolymer.Length == readHomopolymer.Length)
            {
                createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

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
                createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

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

        createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, refMasks, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase, lastMatchedLength, priorRefBase, priorReadBase);

        return realignedVariants;
    }

    @VisibleForTesting
    public static List<UltimaRealignedQualModel> getQualVariants(boolean variantInMergedHomopolymers, final SimpleVariant variant,
            final List<UltimaRealignedQualModel> realignedVariants)
    {
        // TODO: ASK THOMAS. Just returning all variants?
        if(variantInMergedHomopolymers)
        {
            // sandwiched snv/mnv case
            List<UltimaRealignedQualModel> leftInserts = Lists.newArrayList();
            List<UltimaRealignedQualModel> leftDels = Lists.newArrayList();
            int leftIndelBalance = 0;
            List<UltimaRealignedQualModel> rightInserts = Lists.newArrayList();
            List<UltimaRealignedQualModel> rightDels = Lists.newArrayList();
            int rightIndelBalance = 0;
            for(int i = 0; i < realignedVariants.size(); )
            {
                SimpleVariant realignedVariant = realignedVariants.get(i).variant();
                if(realignedVariant.position() < variant.position())
                {
                    leftIndelBalance += realignedVariant.indelLength();
                    if(realignedVariant.isInsert())
                    {
                        leftInserts.add(realignedVariants.get(i));
                    }
                    else
                    {
                        leftDels.add(realignedVariants.get(i));
                    }

                    ++i;
                    continue;
                }

                rightIndelBalance += realignedVariant.indelLength();
                if(realignedVariant.isInsert())
                {
                    rightInserts.add(realignedVariants.get(i));
                }
                else
                {
                    rightDels.add(realignedVariants.get(i));
                }

                ++i;
            }

            List<UltimaRealignedQualModel> qualVariants = Lists.newArrayList();
            if(leftIndelBalance > 0)
            {
                qualVariants.addAll(leftInserts);
            }

            if(leftIndelBalance < 0)
            {
                qualVariants.addAll(leftDels);
            }

            if(rightIndelBalance > 0)
            {
                qualVariants.addAll(rightInserts);
            }

            if(rightIndelBalance < 0)
            {
                qualVariants.addAll(rightDels);
            }

            return qualVariants;
        }

        // non-sandwiched snv/mnv case and indel case
        List<Homopolymer> delHomopolymers;
        List<Homopolymer> insertHomopolymers;
        if(variant.indelLength() > 0)
        {
            delHomopolymers = Lists.newArrayList();
            insertHomopolymers = getHomopolymers(variant.Alt.getBytes(), 1, variant.Alt.length() - 1);
        }
        else if(variant.indelLength() < 0)
        {
            delHomopolymers = getHomopolymers(variant.Ref.getBytes(), 1, variant.Ref.length() - 1);
            insertHomopolymers = Lists.newArrayList();
        }
        else
        {
            delHomopolymers = getHomopolymers(variant.Ref.getBytes(), 0, variant.Ref.length() - 1);
            insertHomopolymers = getHomopolymers(variant.Alt.getBytes(), 0, variant.Alt.length() - 1);
        }

        List<UltimaRealignedQualModel> seqVariants = null;
        List<UltimaRealignedQualModel> leftInserts = Lists.newArrayList();
        List<UltimaRealignedQualModel> leftDels = Lists.newArrayList();
        int leftIndelBalance = 0;
        List<UltimaRealignedQualModel> rightInserts = Lists.newArrayList();
        List<UltimaRealignedQualModel> rightDels = Lists.newArrayList();
        int rightIndelBalance = 0;
        for(int i = 0; i < realignedVariants.size(); )
        {
            UltimaRealignedQualModel realignedVariant = realignedVariants.get(i);
            if(seqVariants == null && i + delHomopolymers.size() + insertHomopolymers.size() - 1 < realignedVariants.size())
            {
                seqVariants = Lists.newArrayList();
                int delIndex = 0;
                int insertIndex = 0;
                while(true)
                {
                    UltimaRealignedQualModel currentVariant = realignedVariants.get(i + delIndex + insertIndex);
                    if(currentVariant.variant().isInsert())
                    {
                        if(insertIndex == insertHomopolymers.size())
                        {
                            seqVariants = null;
                            break;
                        }

                        Homopolymer currentInsert = insertHomopolymers.get(insertIndex);
                        if(currentVariant.variant().indelLengthAbs() == currentInsert.Length
                                && currentVariant.variant().Alt.charAt(1) == (char) currentInsert.Base)
                        {
                            seqVariants.add(currentVariant);
                            ++insertIndex;
                        }
                        else
                        {
                            seqVariants = null;
                            break;
                        }
                    }
                    else
                    {
                        if(delIndex == delHomopolymers.size())
                        {
                            seqVariants = null;
                            break;
                        }

                        Homopolymer currentDel = delHomopolymers.get(delIndex);
                        if(currentVariant.variant().indelLengthAbs() == currentDel.Length
                                && currentVariant.variant().Ref.charAt(1) == (char) currentDel.Base)
                        {
                            seqVariants.add(currentVariant);
                            ++delIndex;
                        }
                        else
                        {
                            seqVariants = null;
                            break;
                        }
                    }

                    if(insertIndex == insertHomopolymers.size() && delIndex == delHomopolymers.size())
                    {
                        break;
                    }
                }

                if(seqVariants != null)
                {
                    i += seqVariants.size();
                    continue;
                }
            }

            if(seqVariants == null)
            {
                leftIndelBalance += realignedVariant.variant().indelLength();
                if(realignedVariant.variant().isInsert())
                {
                    leftInserts.add(realignedVariant);
                }
                else
                {
                    leftDels.add(realignedVariant);
                }

                ++i;
                continue;
            }

            rightIndelBalance += realignedVariant.variant().indelLength();
            if(realignedVariant.variant().isInsert())
            {
                rightInserts.add(realignedVariant);
            }
            else
            {
                rightDels.add(realignedVariant);
            }

            ++i;
        }

        if(seqVariants == null)
        {
            return realignedVariants;
        }

        List<UltimaRealignedQualModel> qualVariants = Lists.newArrayList();
        if(leftIndelBalance > 0)
        {
            qualVariants.addAll(leftInserts);
        }

        if(leftIndelBalance < 0)
        {
            qualVariants.addAll(leftDels);
        }

        qualVariants.addAll(seqVariants);

        if(rightIndelBalance > 0)
        {
            qualVariants.addAll(rightInserts);
        }

        if(rightIndelBalance < 0)
        {
            qualVariants.addAll(rightDels);
        }

        return qualVariants;
    }
}
