package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.internal.Sets;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.pcollections.PVector;
import org.pcollections.TreePVector;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// TODO: LATER Move to a better location.
// TODO: Guard calling this against config
// TODO: clean up unneeded functions.
// TODO: LATER comprehensive unit tests.
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

    private static int getVariantRefIndex(final VariantReadContext readContext)
    {
        // TODO: LATER This is repetitive?
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;

        int refIndex = 0;
        int readIndex = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(readIndex == varCoreIndex)
            {
                break;
            }

            if(op.consumesReferenceBases())
            {
                ++refIndex;
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }
        }

        return refIndex;
    }

    private static boolean isSandwichedNonIndel(final VariantReadContext readContext, int varRefIndex)
    {
        if(readContext.variant().isIndel())
        {
            return false;
        }

        int varLength = readContext.variant().Ref.length();

        // TODO: Check that it isn't flanked by del or ins?
        // TODO: LATER START this is repeat code
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;

        int cigarIndex = 0;
        int readIndex = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(readIndex == varCoreIndex)
            {
                break;
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }

            ++cigarIndex;
        }

        final List<CigarOperator> flankedVariantCigarOperators = coreCigarOps
                .subList(cigarIndex - 1, (cigarIndex - 1) + varLength + 2);

        assert flankedVariantCigarOperators.size() == varLength + 2;

        Set<CigarOperator> collapsedFlankedVariantCigarOps = flankedVariantCigarOperators
                .stream()
                .collect(Collectors.toCollection(() -> EnumSet.noneOf(CigarOperator.class)));

        assert collapsedFlankedVariantCigarOps.size() == 1;
        assert collapsedFlankedVariantCigarOps.stream().findFirst().orElse(null) == CigarOperator.M;
        // TODO: LATER END

        String flankedRef = new String(readContext.RefBases, varRefIndex - 1, varLength + 2);
        String flankedRead = new String(readContext.ReadBases, readContext.VarIndex - 1, varLength + 2);

        assert flankedRef.charAt(0) == flankedRead.charAt(0);
        assert flankedRef.charAt(flankedRef.length() - 1) == flankedRead.charAt(flankedRead.length() - 1);

        if(flankedRef.charAt(0) != flankedRef.charAt(flankedRef.length() - 1))
        {
            return false;
        }

        char flankBase = flankedRef.charAt(0);
        String ref = readContext.variant().Ref;
        String alt = readContext.variant().Alt;
        String repeatedFlankBase = String.valueOf(flankBase).repeat(varLength);

        return ref.equals(repeatedFlankBase) || alt.equals(repeatedFlankBase);
    }

    public static List<CigarOperator> getCoreCigarOps(final VariantReadContext readContext, final List<CigarOperator> cigarElements)
    {
        int cigarIndex = 0;
        int readIndex = 0;
        List<CigarOperator> coreCigarOps = Lists.newArrayList();
        while(readIndex <= readContext.CoreIndexEnd)
        {
            CigarOperator op = cigarElements.get(cigarIndex++);
            boolean inCore = readIndex >= readContext.CoreIndexStart;
            if(inCore)
            {
                coreCigarOps.add(op);
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }
        }

        return coreCigarOps;
    }

    // TODO: LATER Use library method.
    public static List<CigarOperator> cigarStringToOps(final String readCigar)
    {
        List<CigarOperator> cigarOps= Lists.newArrayList();
        int count = 0;
        for(int i = 0; i < readCigar.length(); ++i)
        {
            char c = readCigar.charAt(i);
            if('0' <= c && c <= '9')
            {
                count = 10*count + (int)c - (int)'0';
                continue;
            }

            if(count == 0)
            {
                throw new RuntimeException("TODO: Count should not be zero");
            }

            CigarOperator op;
            switch(c)
            {
                case 'M':
                    op = CigarOperator.M;
                    break;
                case 'D':
                    op = CigarOperator.D;
                    break;
                default:
                    throw new RuntimeException("TODO: Cigar operator not recognized");
            }

            for(; count > 0; --count)
            {
                cigarOps.add(op);
            }
        }

        if(count > 0)
        {
            throw new RuntimeException("TODO: Count should be zero");
        }

        return cigarOps;
    }

    public static List<UltimaQualModel> buildUltimaQualModels(final VariantReadContext readContext, final UltimaQualCalculator ultimaQualCalculator)
    {
        final SimpleVariant variant = readContext.variant();

        // TODO: Test inserts
        assert variant.indelLength() <= 0;

        // TODO: consider mnv variants.
        assert !variant.isMNV();

        // TODO: Consider sandwiched SNV/MNV.
        assert !isSandwichedNonIndel(readContext, getVariantRefIndex(readContext));

        // TODO: Mask sandwich snvs/mvns.

        List<Homopolymer> refHomopolymers = getHomopolymers(readContext.refBases());
        List<Homopolymer> readHomopolymers = getHomopolymers(readContext.coreStr());

        // TODO: extract method?
        // construct realigned indel variants
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

        // TODO: extract method?
        // collect variants that contribute to final qual score
        if(isSandwichedNonIndel(readContext, getVariantRefIndex(readContext)))
        {
            // TODO: Consider sandwiched SNV/MNV.
            assert false;
            return null;
        }

        List<SimpleVariant> qualVariants = Lists.newArrayList();
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
                qualVariants = realignedVariants;
            }
            else
            {
                // TODO: QUESTION q_seq in indel context?
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
            }
        }
        else
        {
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
                qualVariants = realignedVariants;
            }
            else
            {
                // TODO: QUESTION q_seq in MNV context?
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
            }
        }

        // turn variants into a qual model
        return qualVariants.stream().map(x -> ultimaQualCalculator.buildContext(x)).collect(Collectors.toList());
    }
}
