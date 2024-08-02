package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.Nullable;

// TODO: LATER Move to a better location.
// TODO: clean up unneeded functions.
// TODO: LATER comprehensive unit tests.
// TODO: Test on actual sample.
// TODO: LATER performance testing.
public class UltimaRealignedQualModelBuilder
{
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
            return String.valueOf((char) Base) + "x" + Length;
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
        List<SimpleVariant> realignedVariants = getRealignedVariants(readContext, mergedHomopolymers.RefHomopolymers, mergedHomopolymers.ReadHomopolymers);
        List<SimpleVariant> qualVariants = getQualVariants(mergedHomopolymers.variantInMergedHomopolymers(), readContext, realignedVariants);
        return qualVariants.stream().map(x -> ultimaQualCalculator.buildContext(x)).collect(Collectors.toList());
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
            if(homopolymer.Base > prev_homopolymer.Base)
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
                    // TODO: simplify this
                    // look forward in ref
                    // TODO: Test case around this.
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
                    // TODO: simplify this
                    // look forward in ref
                    // TODO: Test case around this.
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

    private static void createRealignedVariants(final List<SimpleVariant> realignedVariants, final VariantReadContext readContext,
            final StringBuilder delBases, final StringBuilder insBases, final int lastMatchedRefPos, final byte lastMatchedBase)
    {
        SimpleVariant variant = readContext.variant();
        if(delBases.length() > 0)
        {
            if(lastMatchedRefPos == -1)
            {
                int variantPos = readContext.CorePositionStart - 1;
                String ref = String.valueOf((char) readContext.RefBaseBeforeCore) + delBases.toString();
                String alt = String.valueOf((char) readContext.ReadBases[readContext.CoreIndexStart - 1]);
                realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));
            }
            else
            {
                String ref = String.valueOf((char) lastMatchedBase) + delBases.toString();
                String alt = String.valueOf((char) lastMatchedBase);
                realignedVariants.add(new SimpleVariant(variant.Chromosome, lastMatchedRefPos, ref, alt));
            }
        }

        if(insBases.length() > 0)
        {
            if(lastMatchedRefPos == -1)
            {
                int variantPos = readContext.CorePositionStart - 1;
                String ref = String.valueOf((char) readContext.RefBaseBeforeCore);
                String alt = String.valueOf((char) readContext.ReadBases[readContext.CoreIndexStart - 1]) + insBases.toString();
                realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));
            }
            else
            {
                String ref = String.valueOf((char) lastMatchedBase);
                String alt = String.valueOf((char) lastMatchedBase) + insBases.toString();
                realignedVariants.add(new SimpleVariant(variant.Chromosome, lastMatchedRefPos, ref, alt));
            }
        }
    }

    private static List<SimpleVariant> getRealignedVariants(final VariantReadContext readContext, final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers)
    {
        SimpleVariant variant = readContext.variant();
        List<SimpleVariant> realignedVariants = Lists.newArrayList();

        byte lastMatchedBase = 0;
        int lastMatchedRefPos = -1;
        int refIndex = 0;
        int readIndex = 0;
        int refPos = readContext.CorePositionStart;
        StringBuilder delBases = new StringBuilder();
        StringBuilder insBases = new StringBuilder();
        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolymer.Base && refHomopolymer.Length == readHomopolymer.Length)
            {
                createRealignedVariants(realignedVariants, readContext, delBases, insBases, lastMatchedRefPos, lastMatchedBase);

                lastMatchedBase = refHomopolymer.Base;
                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                lastMatchedRefPos = refPos - 1;
                delBases = new StringBuilder();
                insBases = new StringBuilder();
                continue;
            }

            if(refHomopolymer.Base == readHomopolymer.Base)
            {
                createRealignedVariants(realignedVariants, readContext, delBases, insBases, lastMatchedRefPos, lastMatchedBase);

                lastMatchedBase = refHomopolymer.Base;
                ++refIndex;
                ++readIndex;
                lastMatchedRefPos = refPos + min(refHomopolymer.Length, readHomopolymer.Length) - 1;
                refPos += refHomopolymer.Length;
                delBases = new StringBuilder();
                insBases = new StringBuilder();

                if(refHomopolymer.Length < readHomopolymer.Length)
                {
                    insBases.append(String.valueOf((char) refHomopolymer.Base).repeat(readHomopolymer.Length - refHomopolymer.Length));
                }
                else
                {
                    delBases.append(String.valueOf((char) refHomopolymer.Base).repeat(refHomopolymer.Length - readHomopolymer.Length));
                }

                continue;
            }

            int refHomopolymersLeft = refHomopolymers.size() - refIndex - 1;
            int readHomopolymersLeft = readHomopolymers.size() - readIndex - 1;
            if(refHomopolymersLeft >= readHomopolymersLeft)
            {
                delBases.append(String.valueOf((char) refHomopolymer.Base).repeat(refHomopolymer.Length));
                ++refIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            insBases.append(String.valueOf((char) readHomopolymer.Base).repeat(readHomopolymer.Length));
            ++readIndex;
        }

        while(refIndex < refHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            delBases.append(String.valueOf((char) refHomopolymer.Base).repeat(refHomopolymer.Length));
            ++refIndex;
        }

        while(readIndex < readHomopolymers.size())
        {
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);
            insBases.append(String.valueOf((char) readHomopolymer.Base).repeat(readHomopolymer.Length));
            ++readIndex;
        }

        createRealignedVariants(realignedVariants, readContext, delBases, insBases, lastMatchedRefPos, lastMatchedBase);

        return realignedVariants;
    }

    private static List<SimpleVariant> getQualVariants(boolean variantInMergedHomopolymers, final VariantReadContext readContext, final List<SimpleVariant> realignedVariants)
    {
        SimpleVariant variant = readContext.variant();

        // TODO: remove this repetition.
        if(variantInMergedHomopolymers)
        {
            // sandwiched snv/mnv case
            // TODO: LATER lots of repetition in terms of left/right indels.
            List<SimpleVariant> leftInserts = Lists.newArrayList();
            List<SimpleVariant> leftDels = Lists.newArrayList();
            int leftIndelBalance = 0;
            List<SimpleVariant> rightInserts = Lists.newArrayList();
            List<SimpleVariant> rightDels = Lists.newArrayList();
            int rightIndelBalance = 0;
            for(int i = 0; i < realignedVariants.size();)
            {
                SimpleVariant realignedVariant = realignedVariants.get(i);
                if(realignedVariant.position() < variant.position())
                {
                    leftIndelBalance += realignedVariant.indelLength();
                    if(realignedVariant.isInsert())
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

                rightIndelBalance += realignedVariant.indelLength();
                if(realignedVariant.isInsert())
                {
                    rightInserts.add(realignedVariant);
                }
                else
                {
                    rightDels.add(realignedVariant);
                }

                ++i;
            }

            List<SimpleVariant> qualVariants = Lists.newArrayList();
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

        // TODO: remove this repetition.
        if(variant.isIndel())
        {
            int indelLength = variant.indelLength();
            char indelBase = indelLength > 0 ? variant.Alt.charAt(1) : variant.Ref.charAt(1);

            int variantIndex = -1;
            List<SimpleVariant> leftInserts = Lists.newArrayList();
            List<SimpleVariant> leftDels = Lists.newArrayList();
            int leftIndelBalance = 0;
            List<SimpleVariant> rightInserts = Lists.newArrayList();
            List<SimpleVariant> rightDels = Lists.newArrayList();
            int rightIndelBalance = 0;
            for(int i = 0; i < realignedVariants.size(); ++i)
            {
                SimpleVariant realignedVariant = realignedVariants.get(i);
                int realignedIndelLength = realignedVariant.indelLength();
                char realignedIndelBase = realignedIndelLength > 0 ? realignedVariant.Alt.charAt(1) : realignedVariant.Ref.charAt(1);
                if(variantIndex == -1
                        && indelLength == realignedIndelLength
                        && indelBase == realignedIndelBase)
                {
                    variantIndex = i;
                    continue;
                }

                if(variantIndex == -1)
                {
                    leftIndelBalance += realignedVariant.indelLength();
                    if(realignedVariant.isInsert())
                    {
                        leftInserts.add(realignedVariant);
                    }
                    else
                    {
                        leftDels.add(realignedVariant);
                    }

                    continue;
                }

                rightIndelBalance += realignedVariant.indelLength();
                if(realignedVariant.isInsert())
                {
                    rightInserts.add(realignedVariant);
                }
                else
                {
                    rightDels.add(realignedVariant);
                }
            }

            if(variantIndex == -1)
            {
                return realignedVariants;
            }

            List<SimpleVariant> qualVariants = Lists.newArrayList();
            qualVariants.add(realignedVariants.get(variantIndex));
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

        // non-sandwiched snv/mnv case
        List<Homopolymer> delHomopolymers = getHomopolymers(variant.Ref.getBytes(), 0, variant.Ref.length() - 1);
        List<Homopolymer> insertHomopolymers = getHomopolymers(variant.Alt.getBytes(), 0, variant.Alt.length() - 1);

        // TODO: LATER lots of repetition in terms of left/right indels.
        List<SimpleVariant> seqVariants = null;
        List<SimpleVariant> leftInserts = Lists.newArrayList();
        List<SimpleVariant> leftDels = Lists.newArrayList();
        int leftIndelBalance = 0;
        List<SimpleVariant> rightInserts = Lists.newArrayList();
        List<SimpleVariant> rightDels = Lists.newArrayList();
        int rightIndelBalance = 0;
        for(int i = 0; i < realignedVariants.size();)
        {
            SimpleVariant realignedVariant = realignedVariants.get(i);
            if(seqVariants == null && i + delHomopolymers.size() + insertHomopolymers.size() - 1 < realignedVariants.size())
            {
                seqVariants = Lists.newArrayList();
                int delIndex = 0;
                int insertIndex = 0;
                while(true)
                {
                    SimpleVariant currentVariant = realignedVariants.get(i + delIndex + insertIndex);
                    if(currentVariant.isInsert())
                    {
                        if(insertIndex == insertHomopolymers.size())
                        {
                            seqVariants = null;
                            break;
                        }

                        Homopolymer currentInsert = insertHomopolymers.get(insertIndex);
                        if(currentVariant.indelLengthAbs() == currentInsert.Length && currentVariant.Alt.charAt(1) == (char) currentInsert.Base)
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
                        if(currentVariant.indelLengthAbs() == currentDel.Length && currentVariant.Ref.charAt(1) == (char) currentDel.Base)
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
                leftIndelBalance += realignedVariant.indelLength();
                if(realignedVariant.isInsert())
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

            rightIndelBalance += realignedVariant.indelLength();
            if(realignedVariant.isInsert())
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

        List<SimpleVariant> qualVariants = Lists.newArrayList();
        qualVariants.addAll(seqVariants);
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
}
