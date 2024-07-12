package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.String.format;

import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// TODO: Integrate into classes where appropriate
// TODO: Guard calling this against config
// TODO: clean up unneeded functions.
public class UltimaLocalRealigner
{
    private static class Homopolymer
    {
        // TODO: Just use byte instead of char.
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

    private static class HomopolymerPair
    {
        @Nullable
        public final Homopolymer RefHomopolymer;
        @Nullable
        public final Homopolymer ReadHomopolymer;

        public HomopolymerPair(@Nullable final Homopolymer refHomopolymer, @Nullable final Homopolymer readHomopolymer)
        {
            RefHomopolymer = refHomopolymer;
            ReadHomopolymer = readHomopolymer;
        }

        public char base()
        {
            return RefHomopolymer == null ? ReadHomopolymer.Base : RefHomopolymer.Base;
        }

        public int indelLength()
        {
            if(RefHomopolymer != null && ReadHomopolymer != null)
            {
                return ReadHomopolymer.Length - RefHomopolymer.Length;
            }

            return RefHomopolymer == null ? ReadHomopolymer.Length : -RefHomopolymer.Length;
        }

        @Override
        public String toString()
        {
            String refString = RefHomopolymer == null ? "null" : RefHomopolymer.toString();
            String readString = ReadHomopolymer == null ? "null" : ReadHomopolymer.toString();
            return format("ref(%s) read(%s)", refString, readString);
        }
    }

    // TODO: Avoid strings?
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

    // TODO: Does this ever return false? Test this.
    public static boolean isVariantExplained(final VariantReadContext readContext)
    {
        // TODO: consider indel variants.
        final SimpleVariant variant = readContext.variant();

        // TODO: Test longer dels
        // TODO: Test inserts
        assert variant.indelLength() <= 0;

        // TODO: consider mnv variants.
        assert !variant.isMNV();

        // TODO: Consider sandwiched SNV/MNV.
        assert !isSandwichedNonIndel(readContext, getVariantRefIndex(readContext));

        // TODO: Extend so that read base count matches core length modulo indel size, plus one base padding.
        // TODO: Extend read bases and core to cover full homopolymers.
        // TODO: Mask sandwich snvs/mvns.

        List<Homopolymer> refCoreHomopolymers = getHomopolymers(readContext.refBases());
        List<Homopolymer> readCoreHomopolymers = getHomopolymers(readContext.coreStr());

        // Construct read cigar in core
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);
//        List<CigarElement> coreCigarElements = collapseCigarOps(coreCigarOps);

        // TODO: Check core length.
        int actualCoreLength = readContext.coreLength();
        int expectedCoreLength = coreCigarOps.stream().mapToInt(op -> op.consumesReadBases() ? 1 : 0).sum();
        assert actualCoreLength == expectedCoreLength;

        // TODO: Check ref bases length.
        int actualRefBasesLength = readContext.RefBases.length;
        int expectedRefBasesLength = coreCigarOps.stream().mapToInt(op -> op.consumesReferenceBases() ? 1 : 0).sum();
        assert actualRefBasesLength == expectedRefBasesLength;

        // pair homopolymers between ref and read
        List<HomopolymerPair> homopolymerPairs = pairHomopolymers(refCoreHomopolymers, readCoreHomopolymers);

        // TODO: Understand ultima qual calculator. Compute qual. Only if variant still exists after re-alignment.

        return isVariantExplainedHelper(readContext, homopolymerPairs);
    }

    // TODO: come up with a better way
    // TODO: Really need to come up with a better method.
    private static boolean isVariantExplainedHelper(final VariantReadContext readContext, final List<HomopolymerPair> homopolymerPairs)
    {
        SimpleVariant variant = readContext.variant();

        // TODO:
        assert variant.indelLength() <= 0;
        assert !variant.isMNV();

        // now check if homopolymers indels explain the variant.
        if(variant.isSNV())
        {
            // get index of homopolymer pair associated to variant.
            int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart + 1;
            int varHomopolymerPairIndex = -1;

            int readBasesConsumed = 0;
            while(readBasesConsumed < varCoreIndex)
            {
                Homopolymer readHomopolymer = homopolymerPairs.get(++varHomopolymerPairIndex).ReadHomopolymer;
                if(readHomopolymer != null)
                {
                    readBasesConsumed += readHomopolymer.Length;
                }
            }

            HomopolymerPair varHomopolymerPair = homopolymerPairs.get(varHomopolymerPairIndex);
            HomopolymerPair beforeHomopolymerPair = varHomopolymerPairIndex == 0 ? null : homopolymerPairs.get(varHomopolymerPairIndex - 1);
            HomopolymerPair afterHomopolymerPair = varHomopolymerPairIndex == homopolymerPairs.size() - 1 ? null : homopolymerPairs.get(varHomopolymerPairIndex + 1);

            if (varHomopolymerPair.base() != variant.Alt.charAt(0) || varHomopolymerPair.indelLength() != 1)
            {
                return false;
            }

            // check var pair and before pair
            if(beforeHomopolymerPair != null
                    && beforeHomopolymerPair.base() == variant.Ref.charAt(0)
                    && beforeHomopolymerPair.indelLength() == -1)
            {
                return true;
            }

            // check var pair and after pair
            if(afterHomopolymerPair != null
                    && afterHomopolymerPair.base() == variant.Ref.charAt(0)
                    && afterHomopolymerPair.indelLength() == -1)
            {
                return true;
            }

            return false;
        }
        else if(variant.indelLength() <= -1)
        {
            String deletedBases = variant.Ref.substring(1);
            List<Homopolymer> deletedAsHomopolymers = getHomopolymers(deletedBases);

            char deletedBase = deletedBases.charAt(0);
            int matchingPairCount = 0;
            for(int i = 0; i < homopolymerPairs.size(); ++i)
            {
                boolean matches = true;
                for(int j = 0; j < deletedAsHomopolymers.size(); ++j)
                {
                    Homopolymer deleted = deletedAsHomopolymers.get(j);
                    HomopolymerPair pair = homopolymerPairs.get(i + j);
                    if(deleted.Base != pair.base())
                    {
                        matches = false;
                        break;
                    }

                    if(-deleted.Length != pair.indelLength())
                    {
                        matches = false;
                        break;
                    }

                    if(j > 0 && j < deletedAsHomopolymers.size() - 1 && pair.ReadHomopolymer != null)
                    {
                        matches = false;
                        break;
                    }
                }

                if(matches)
                {
                    ++matchingPairCount;
                }
            }

            // TODO: not more than one match.
            assert matchingPairCount <= 1;

            return matchingPairCount == 1;
        }
        else
        {
            // TODO:
            throw new RuntimeException("TODO");
        }
    }

    private static int getVariantRefIndex(final VariantReadContext readContext)
    {
        // TODO: This is repetitive?
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        String coreRefBases = readContext.refBases();
        String coreReadBases = readContext.coreStr();
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
        // TODO: START this is repeat code
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        String coreRefBases = readContext.refBases();
        String coreReadBases = readContext.coreStr();
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
        // TODO: END

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


    // TODO: Use dynamic programming algorithm.
    private static List<HomopolymerPair> pairHomopolymers(final List<Homopolymer> refCoreHomopolymers, final List<Homopolymer> readCoreHomopolymers)
    {
        List<HomopolymerPair> homopolymerPairs = Lists.newArrayList();
        int refIndex = 0;
        int readIndex = 0;
        while(refIndex < refCoreHomopolymers.size() && readIndex < readCoreHomopolymers.size())
        {
            Homopolymer refHomopolymer = refCoreHomopolymers.get(refIndex);
            Homopolymer readHomopolymer = readCoreHomopolymers.get(readIndex);

            if(refHomopolymer.Base == readHomopolymer.Base)
            {
                homopolymerPairs.add(new HomopolymerPair(refHomopolymer, readHomopolymer));
                ++refIndex;
                ++readIndex;
                continue;
            }

            int refRemaining = refCoreHomopolymers.size() - refIndex - 1;
            int readRemaining = readCoreHomopolymers.size() - readIndex - 1;
            if(refRemaining == readRemaining)
            {
                // TODO:
                throw new RuntimeException("TODO");
            }

            if(refRemaining < readRemaining)
            {
                homopolymerPairs.add(new HomopolymerPair(null, readHomopolymer));
                ++readIndex;
                continue;
            }

            homopolymerPairs.add(new HomopolymerPair(refHomopolymer, null));
            ++refIndex;
        }

        while(refIndex < refCoreHomopolymers.size())
        {
            Homopolymer refHomopolymer = refCoreHomopolymers.get(refIndex);
            homopolymerPairs.add(new HomopolymerPair(refHomopolymer, null));
            ++refIndex;
        }

        while(readIndex < readCoreHomopolymers.size())
        {
            Homopolymer readHomopolymer = readCoreHomopolymers.get(readIndex);
            homopolymerPairs.add(new HomopolymerPair(null, readHomopolymer));
            ++readIndex;
        }

        return homopolymerPairs;
    }

    private static List<CigarElement> collapseCigarOps(final List<CigarOperator> coreCigarOps)
    {
        List<CigarElement> cigarElements = Lists.newArrayList();
        CigarOperator currentOp = null;
        int length = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(currentOp == null)
            {
                currentOp = op;
                length = 1;
                continue;
            }

            if(currentOp == op)
            {
                ++length;
                continue;
            }

            cigarElements.add(new CigarElement(length, currentOp));
            currentOp = op;
            length = 1;
        }

        if(length >= 1)
        {
            cigarElements.add(new CigarElement(length, currentOp));
        }

        return cigarElements;
    }

    private static List<CigarOperator> getCoreCigarOps(final VariantReadContext readContext, final List<CigarOperator> cigarElements)
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

    // TODO: Use library method.
    private static List<CigarOperator> cigarStringToOps(final String readCigar)
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
}
