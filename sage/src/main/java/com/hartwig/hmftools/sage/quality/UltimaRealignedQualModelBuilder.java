package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.CYCLE_BASES;

import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

// TODO: clean up unneeded functions.
// TODO: run "Reformat code"
// TODO: comprehensive unit tests.
public class UltimaRealignedQualModelBuilder
{
    private static final HashMap<Byte, Integer> CYCLE_BASE_INDEX;
    static
    {
        CYCLE_BASE_INDEX = Maps.newHashMap();
        for(int i = 0; i < CYCLE_BASES.length; i++)
        {
            CYCLE_BASE_INDEX.put(CYCLE_BASES[i], i);
        }
    }

    // TODO: This is duct tape for now.
    public static class UltimaRealignedQualModel extends UltimaQualModel
    {
        private final SimpleVariant mVariant;
        private final UltimaQualModel mBaseQualModel;
        private final int mVarReadIndexOffset;

        // TODO: Just need offset, because core will always match right?
        public UltimaRealignedQualModel(final SimpleVariant variant, final UltimaQualModel baseQualModel, int varReadIndexOffset)
        {
            super(baseQualModel.type());

            mVariant = variant;
            mBaseQualModel = baseQualModel;
            mVarReadIndexOffset = varReadIndexOffset;
        }

        @VisibleForTesting
        public UltimaRealignedQualModel(final SimpleVariant variant)
        {
            // this is just a placeholder for testing
            super(UltimaModelType.HOMOPOLYMER_ADJUSTMENT);

            mVariant = variant;
            mBaseQualModel = null;
            mVarReadIndexOffset = -1;
        }

        @Override
        public byte calculateQual(final SAMRecord record, final int varReadIndex)
        {
            return mBaseQualModel.calculateQual(record, varReadIndex + mVarReadIndexOffset);
        }

        public SimpleVariant variant() { return mVariant; }
    }

    @VisibleForTesting
    public static class Homopolymer
    {
        public final byte Base;
        public final int Length;

        public Homopolymer(final byte base, final int length)
        {
            Base = base;
            Length = length;
        }

        @Override
        public String toString()
        {
            return String.valueOf(Length) + "x" + (char) Base;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof Homopolymer))
            {
                return false;
            }
            final Homopolymer that = (Homopolymer) o;
            return Base == that.Base && Length == that.Length;
        }

        @Override
        public int hashCode()
        {
            return (int) Base + 31*Length;
        }
    }

    public static List<UltimaQualModel> buildRealignedUltimaQualModels(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator)
    {
        List<Homopolymer> refHomopolymers = getHomopolymers(readContext.RefBases, 0, readContext.RefBases.length - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers);
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(readContext, ultimaQualCalculator, mergedHomopolymers.RefHomopolymers, mergedHomopolymers.ReadHomopolymers);

        // TODO: Move this into a unit test.
        for(int i = 1; i < realignedVariants.size(); i++)
        {
            if(realignedVariants.get(i).variant().Position < realignedVariants.get(i - 1).variant().Position)
            {
                throw new IllegalStateException("Realigned variants are out of order.");
            }
        }

        // TODO: This has been duct taped.
        List<UltimaQualModel> realignedQualModels =
                getQualVariants(mergedHomopolymers.variantInMergedHomopolymers(), readContext.variant(), realignedVariants)
                        .stream()
                        .map(x -> (UltimaQualModel) x)
                        .collect(Collectors.toList());

        if(realignedQualModels.isEmpty() && !mergedHomopolymers.variantInMergedHomopolymers())
        {
            throw new IllegalStateException("Variant({}) is expected to have realigned ultima variants, but none have been found.");
        }

        return realignedQualModels;
    }

    @VisibleForTesting
    public static List<Homopolymer> getHomopolymers(final byte[] bases, int startIndex, int endIndex)
    {
        List<Homopolymer> homopolymers = Lists.newArrayList();
        if(endIndex < startIndex)
        {
            return homopolymers;
        }

        byte currentBase = bases[startIndex];
        int currentLength = 1;
        for(int i = startIndex + 1; i <= endIndex; i++)
        {
            byte base = bases[i];
            if(base == currentBase)
            {
                currentLength++;
                continue;
            }

            homopolymers.add(new Homopolymer(currentBase, currentLength));
            currentBase = base;
            currentLength = 1;
        }

        homopolymers.add(new Homopolymer(currentBase, currentLength));
        return homopolymers;
    }

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

    @VisibleForTesting
    public static class MergedHomopolymers
    {
        public final List<Homopolymer> RefHomopolymers;
        public final List<Homopolymer> ReadHomopolymers;
        private boolean VariantInMergedHomopolymers;

        private MergedHomopolymers(final List<Homopolymer> refHomopolymers, final List<Homopolymer> readHomopolymers, boolean variantInMergedHomopolymers)
        {
            RefHomopolymers = refHomopolymers;
            ReadHomopolymers = readHomopolymers;
            VariantInMergedHomopolymers = variantInMergedHomopolymers;
        }

        public MergedHomopolymers()
        {
            this(Lists.newArrayList(), Lists.newArrayList(), false);
        }

        public void setVariantInMergedHomopolymers()
        {
            VariantInMergedHomopolymers = true;
        }

        public boolean variantInMergedHomopolymers()
        {
            return VariantInMergedHomopolymers;
        }
    }

    // TODO: Simplify this.
    @VisibleForTesting
    public static MergedHomopolymers mergeSandwichedHomopolymers(@Nullable final VariantReadContext readContext,
            final List<Homopolymer> refHomopolymers0, final List<Homopolymer> readHomopolymers0)
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(refHomopolymers0);
        List<Homopolymer> readHomopolymers = Lists.newArrayList(readHomopolymers0);

        MergedHomopolymers mergedHomopolymers = new MergedHomopolymers();
        int refIndex = 0;
        int readIndex = 0;
        int readBasesConsumed = 0;
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
                        int refCycleCount = cycleCount(refHomopolymers, refIndex);
                        int readCycleCount = cycleCount(readHomopolymers, readIndex);
                        while(true)
                        {
                            if(forwardRefHompolymer.Base == refHomopolymer.Base && refCycleCount != readCycleCount)
                            {
                                int mergedLength = IntStream.range(refIndex, refIndex + extraRefIndex + 1).map(i -> refHomopolymers.get(i).Length).sum();
                                refHomopolymers.set(refIndex, new Homopolymer(refHomopolymer.Base, mergedLength));
                                for(int i = refIndex + 1; i <= refIndex + extraRefIndex; i++)
                                {
                                    refHomopolymers.remove(refIndex + 1);
                                }
                                contracted = true;

                                if(readContext != null && !readContext.variant().isIndel())
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
                            if(forwardReadHompolymer.Base == refHomopolymer.Base && refCycleCount != readCycleCount)
                            {
                                int mergedLength = IntStream.range(readIndex, readIndex + extraReadIndex + 1).map(i -> readHomopolymers.get(i).Length).sum();
                                readHomopolymers.set(readIndex, new Homopolymer(refHomopolymer.Base, mergedLength));
                                for(int i = readIndex + 1; i <= readIndex + extraReadIndex; i++)
                                {
                                    readHomopolymers.remove(readIndex + 1);
                                }
                                contracted = true;

                                if(readContext != null && !readContext.variant().isIndel())
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
                        readIndex++;
                        readBasesConsumed += readHomopolymer.Length;
                    }
                }
                else
                {
                    mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                    mergedHomopolymers.ReadHomopolymers.add(readHomopolymer);
                    refIndex++;
                    readIndex++;
                    readBasesConsumed += readHomopolymer.Length;
                }
            }
            else if(readRemaining <= refRemaining)
            {
                // novel delete
                mergedHomopolymers.RefHomopolymers.add(refHomopolymer);
                refIndex++;
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

    private static void createRealignedVariants(final List<UltimaRealignedQualModel> realignedVariants,
            final UltimaQualCalculator ultimaQualCalculator, final VariantReadContext readContext, final List<Homopolymer> delHomopolymers,
            final List<Homopolymer> insHomopolymers, final int lastMatchedReadCoreIndex, final int lastMatchedRefPos, final byte lastMatchedBase)
    {
        if(delHomopolymers.isEmpty() && insHomopolymers.isEmpty())
        {
            return;
        }

        String chromosome = readContext.variant().chromosome();
        if(lastMatchedRefPos == -1 || lastMatchedReadCoreIndex == -1)
        {
            throw new IllegalArgumentException("Looks as though the first core homopolymers of the read and ref do not have the same base.");
        }

        int origReadCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;
        int varReadIndexOffset = lastMatchedReadCoreIndex - origReadCoreIndex;
        for(int i = 0; i < delHomopolymers.size(); i++)
        {
            Homopolymer delHomopolymer = delHomopolymers.get(i);
            String delBasesString = String.valueOf((char) delHomopolymer.Base).repeat(delHomopolymer.Length);
            String ref = String.valueOf((char) lastMatchedBase) + delBasesString;
            String alt = String.valueOf((char) lastMatchedBase);

            SimpleVariant variant = new SimpleVariant(chromosome, lastMatchedRefPos, ref, alt);
            UltimaQualModel baseQualModel = ultimaQualCalculator.buildContext(variant);
            realignedVariants.add(new UltimaRealignedQualModel(variant, baseQualModel, varReadIndexOffset));
        }

        for(int i = 0; i < insHomopolymers.size(); i++)
        {
            Homopolymer insHomopolymer = insHomopolymers.get(i);
            String insBasesString = String.valueOf((char) insHomopolymer.Base).repeat(insHomopolymer.Length);
            String ref = String.valueOf((char) lastMatchedBase);
            String alt = String.valueOf((char) lastMatchedBase) + insBasesString;

            SimpleVariant variant = new SimpleVariant(chromosome, lastMatchedRefPos, ref, alt);
            UltimaQualModel baseQualModel = ultimaQualCalculator.buildContext(variant);
            realignedVariants.add(new UltimaRealignedQualModel(variant, baseQualModel, varReadIndexOffset));
        }
    }

    private static void extendHomopolymers(final List<Homopolymer> homopolymers, byte base, int length)
    {
        if(homopolymers.isEmpty())
        {
            homopolymers.add(new Homopolymer(base, length));
            return;
        }

        Homopolymer lastHomopolymer = homopolymers.get(homopolymers.size() - 1);
        if(lastHomopolymer.Base == base)
        {
            homopolymers.set(homopolymers.size() - 1, new Homopolymer(base, lastHomopolymer.Length + length));
            return;
        }

        homopolymers.add(new Homopolymer(base, length));
    }

    private static List<UltimaRealignedQualModel> getRealignedVariants(final VariantReadContext readContext,
            final UltimaQualCalculator ultimaQualCalculator, final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers)
    {
        List<UltimaRealignedQualModel> realignedVariants = Lists.newArrayList();

        int lastMatchedReadCoreIndex = -1;
        int readCoreIndex = 0;
        byte lastMatchedBase = 0;
        int lastMatchedRefPos = -1;
        int refIndex = 0;
        int readIndex = 0;
        int refPos = readContext.CorePositionStart;
        List<Homopolymer> delHomopolymers = Lists.newArrayList();
        List<Homopolymer> insHomopolymers = Lists.newArrayList();
        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolymer.Base && refHomopolymer.Length == readHomopolymer.Length)
            {
                createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase);

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
                createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase);

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
                    insHomopolymers.add(new Homopolymer(refHomopolymer.Base, readHomopolymer.Length - refHomopolymer.Length));
                }
                else
                {
                    delHomopolymers.add(new Homopolymer(refHomopolymer.Base, refHomopolymer.Length - readHomopolymer.Length));
                }

                continue;
            }

            int refHomopolymersLeft = refHomopolymers.size() - refIndex - 1;
            int readHomopolymersLeft = readHomopolymers.size() - readIndex - 1;
            if(refHomopolymersLeft == readHomopolymersLeft)
            {
                extendHomopolymers(delHomopolymers, refHomopolymer.Base, refHomopolymer.Length);
                extendHomopolymers(insHomopolymers, readHomopolymer.Base, readHomopolymer.Length);
                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                readCoreIndex += readHomopolymer.Length;
                continue;
            }

            if(refHomopolymersLeft > readHomopolymersLeft)
            {
                extendHomopolymers(delHomopolymers, refHomopolymer.Base, refHomopolymer.Length);
                ++refIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            extendHomopolymers(insHomopolymers, readHomopolymer.Base, readHomopolymer.Length);
            ++readIndex;
            readCoreIndex += readHomopolymer.Length;
        }

        while(refIndex < refHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            extendHomopolymers(delHomopolymers, refHomopolymer.Base, refHomopolymer.Length);
            ++refIndex;
        }

        while(readIndex < readHomopolymers.size())
        {
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);
            extendHomopolymers(insHomopolymers, readHomopolymer.Base, readHomopolymer.Length);
            ++readIndex;
        }

        createRealignedVariants(realignedVariants, ultimaQualCalculator, readContext, delHomopolymers, insHomopolymers, lastMatchedReadCoreIndex, lastMatchedRefPos, lastMatchedBase);

        return realignedVariants;
    }

    // TODO: Simplify and remove repetition.
    @VisibleForTesting
    public static List<UltimaRealignedQualModel> getQualVariants(boolean variantInMergedHomopolymers, final SimpleVariant variant, final List<UltimaRealignedQualModel> realignedVariants)
    {
        if(variantInMergedHomopolymers)
        {
            // sandwiched snv/mnv case
            List<UltimaRealignedQualModel> leftInserts = Lists.newArrayList();
            List<UltimaRealignedQualModel> leftDels = Lists.newArrayList();
            int leftIndelBalance = 0;
            List<UltimaRealignedQualModel> rightInserts = Lists.newArrayList();
            List<UltimaRealignedQualModel> rightDels = Lists.newArrayList();
            int rightIndelBalance = 0;
            for(int i = 0; i < realignedVariants.size();)
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
        for(int i = 0; i < realignedVariants.size();)
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
                        if(currentVariant.variant().indelLengthAbs() == currentInsert.Length && currentVariant.variant().Alt.charAt(1) == (char) currentInsert.Base)
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
                        if(currentVariant.variant().indelLengthAbs() == currentDel.Length && currentVariant.variant().Ref.charAt(1) == (char) currentDel.Base)
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
