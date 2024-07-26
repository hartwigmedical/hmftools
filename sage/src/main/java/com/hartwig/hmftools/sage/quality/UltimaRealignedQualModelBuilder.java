package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;

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
        public final char Base;
        public final int Length;

        public Homopolymer(final char base, final int length)
        {
            Base = base;
            Length = length;
        }

        @Override
        public String toString()
        {
            return String.valueOf(Base) + "x" + Length;
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
            return Objects.hash(Base, Length);
        }
    }

    public static List<UltimaQualModel> buildRealignedUltimaQualModels(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator)
    {
        // TODO: QUESTION Consider sandwiched SNV/MNV.
        // TODO: QUESTION Mask sandwich snvs/mvns.

        List<Homopolymer> refHomopolymers = getHomopolymers(readContext.refBases());
        List<Homopolymer> readHomopolymers = getHomopolymers(readContext.coreStr());
        List<SimpleVariant> realignedVariants = getRealignedVariants(readContext, refHomopolymers, readHomopolymers);
        List<SimpleVariant> qualVariants = getQualVariants(readContext, realignedVariants);
        return qualVariants.stream().map(x -> ultimaQualCalculator.buildContext(x)).collect(Collectors.toList());
    }

    private static List<Homopolymer> getHomopolymers(final String bases)
    {
        List<Homopolymer> homopolymers = Lists.newArrayList();
        if(bases.length() == 0)
        {
            return homopolymers;
        }

        char currentBase = bases.charAt(0);
        int currentLength = 1;
        for(int i = 1; i < bases.length(); i++)
        {
            char base = bases.charAt(i);
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

    private static void createRealignedVariants(final List<SimpleVariant> realignedVariants, final VariantReadContext readContext,
            final StringBuilder delBases, final StringBuilder insBases, final int lastMatchedRefPos, final char lastMatchedBase)
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
                String ref = String.valueOf(lastMatchedBase) + delBases.toString();
                String alt = String.valueOf(lastMatchedBase);
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
                String ref = String.valueOf(lastMatchedBase);
                String alt = String.valueOf(lastMatchedBase) + insBases.toString();
                realignedVariants.add(new SimpleVariant(variant.Chromosome, lastMatchedRefPos, ref, alt));
            }
        }
    }

    private static List<SimpleVariant> getRealignedVariants(final VariantReadContext readContext, final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers)
    {
        SimpleVariant variant = readContext.variant();
        List<SimpleVariant> realignedVariants = Lists.newArrayList();

        char lastMatchedBase = (char) 0;
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
                    insBases.append(String.valueOf(refHomopolymer.Base).repeat(readHomopolymer.Length - refHomopolymer.Length));
                }
                else
                {
                    delBases.append(String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length - readHomopolymer.Length));
                }

                continue;
            }

            int refHomopolymersLeft = refHomopolymers.size() - refIndex - 1;
            int readHomopolymersLeft = readHomopolymers.size() - readIndex - 1;
            if(refHomopolymersLeft > readHomopolymersLeft)
            {
                delBases.append(String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length));
                ++refIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            if(refHomopolymersLeft < readHomopolymersLeft)
            {
                insBases.append(String.valueOf(readHomopolymer.Base).repeat(readHomopolymer.Length));
                ++readIndex;
                continue;
            }

            delBases.append(String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length));
            insBases.append(String.valueOf(readHomopolymer.Base).repeat(readHomopolymer.Length));
            ++refIndex;
            ++readIndex;
            refPos += refHomopolymer.Length;
        }

        while(refIndex < refHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            delBases.append(String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length));
            ++refIndex;
        }

        while(readIndex < readHomopolymers.size())
        {
            Homopolymer readHomopolymer = readHomopolymers.get(readIndex);
            insBases.append(String.valueOf(readHomopolymer.Base).repeat(readHomopolymer.Length));
            ++readIndex;
        }

        createRealignedVariants(realignedVariants, readContext, delBases, insBases, lastMatchedRefPos, lastMatchedBase);

        return realignedVariants;
    }

    private static List<SimpleVariant> getQualVariants(final VariantReadContext readContext, final List<SimpleVariant> realignedVariants)
    {
        SimpleVariant variant = readContext.variant();

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

            // TODO: QUESTION q_seq in indel context?
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

        // snv/mnv case
        List<Homopolymer> delHomopolymers = getHomopolymers(variant.Ref);
        List<Homopolymer> insertHomopolymers = getHomopolymers(variant.Alt);

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
                        if(currentVariant.indelLengthAbs() == currentInsert.Length && currentVariant.Alt.charAt(1) == currentInsert.Base)
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
                        if(currentVariant.indelLengthAbs() == currentDel.Length && currentVariant.Ref.charAt(1) == currentDel.Base)
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

        // TODO: QUESTION q_seq in MNV context?
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
