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

    private static List<SimpleVariant> getRealignedVariants(final VariantReadContext readContext, final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers)
    {
        SimpleVariant variant = readContext.variant();
        List<SimpleVariant> realignedVariants = Lists.newArrayList();
        int refIndex = 0;
        int readIndex = 0;
        int refPos = readContext.CorePositionStart;
        while(refIndex < refHomopolymers.size() && readIndex < readHomopolymers.size())
        {
            Homopolymer refHomopolymer = refHomopolymers.get(refIndex);
            Homopolymer readHomopolyer = readHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolyer.Base && refHomopolymer.Length == readHomopolyer.Length)
            {
                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            if(refHomopolymer.Base == readHomopolyer.Base)
            {
                if(refHomopolymer.Length < readHomopolyer.Length)
                {
                    // ins
                    int variantPos = refPos + refHomopolymer.Length - 1;
                    String ref = String.valueOf(refHomopolymer.Base);
                    String alt = ref.repeat(readHomopolyer.Length - refHomopolymer.Length + 1);
                    realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));
                }
                else
                {
                    // del
                    int variantPos = refPos + readHomopolyer.Length - 1;
                    String alt = String.valueOf(refHomopolymer.Base);
                    String ref = alt.repeat(refHomopolymer.Length - readHomopolyer.Length + 1);
                    realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));
                }

                ++refIndex;
                ++readIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            int refHomopolymersLeft = refHomopolymers.size() - refIndex - 1;
            int readHomopolymersLeft = readHomopolymers.size() - readIndex - 1;
            if(refHomopolymersLeft == readHomopolymersLeft)
            {
                // first create novel del from ref homopolymer
                int variantPos = refPos - 1;
                int refBasesIndex = refPos - readContext.CorePositionStart;
                char baseBefore = (char) readContext.RefBases[refBasesIndex - 1];
                String ref = String.valueOf(baseBefore) + String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length);
                String alt = String.valueOf(baseBefore);
                realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));

                refPos += refHomopolymer.Length;

                // now create novel insert from read homopolymer
                variantPos = refPos - 1;
                ref = String.valueOf(refHomopolymer.Base);
                alt = String.valueOf(refHomopolymer.Base) + String.valueOf(readHomopolyer.Base).repeat(readHomopolyer.Length);
                realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));

                ++refIndex;
                ++readIndex;
                continue;
            }

            if(refHomopolymersLeft > readHomopolymersLeft)
            {
                // create novel del from ref homopolymer
                int variantPos = refPos - 1;
                int refBasesIndex = refPos - readContext.CorePositionStart;
                char baseBefore = (char) readContext.RefBases[refBasesIndex - 1];
                String ref = String.valueOf(baseBefore) + String.valueOf(refHomopolymer.Base).repeat(refHomopolymer.Length);
                String alt = String.valueOf(baseBefore);
                realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));

                ++refIndex;
                refPos += refHomopolymer.Length;
                continue;
            }

            // create novel ins from read homopolymer
            int variantPos = refPos - 1;
            int refBasesIndex = refPos - readContext.CorePositionStart;
            char baseBefore = (char) readContext.RefBases[refBasesIndex - 1];
            String ref = String.valueOf(baseBefore);
            String alt = String.valueOf(baseBefore) + String.valueOf(readHomopolyer.Base).repeat(readHomopolyer.Length);
            realignedVariants.add(new SimpleVariant(variant.Chromosome, variantPos, ref, alt));

            ++readIndex;
        }

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
