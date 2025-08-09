package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.CYCLE_BASES;

import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.common.VariantReadContext;

public class UltimaRealignedQualModelsBuilder
{
    // TODO: Access modifiers.
    static final byte INVALID_BASE = -1;

    static final HashMap<Byte, Integer> CYCLE_BASE_INDEX;

    static
    {
        CYCLE_BASE_INDEX = Maps.newHashMap();
        for(int i = 0; i < CYCLE_BASES.length; i++)
        {
            CYCLE_BASE_INDEX.put(CYCLE_BASES[i], i);
        }
    }

    public static class Homopolymer
    {
        public final byte Base;
        public final int Length;

        public Homopolymer(byte base, int length)
        {
            Base = base;
            Length = length;
        }

        public String expand()
        {
            return String.valueOf((char) Base).repeat(Length);
        }

        @Override
        public String toString()
        {
            return Length + "x" + (char) Base;
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
            return (int) Base + 31 * Length;
        }
    }

    public static List<Homopolymer> getHomopolymers(final byte[] bases, int startIndex, int endIndex)
    {
        List<Homopolymer> homopolymers = Lists.newArrayList();
        if(startIndex < 0 || endIndex >= bases.length || endIndex < startIndex)
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

    public static class RefMask
    {
        public final int PosStart;
        public final int PosEnd;
        public final byte BaseMask;

        public RefMask(int posStart, int posEnd, byte baseMask)
        {
            PosStart = posStart;
            PosEnd = posEnd;
            BaseMask = baseMask;
        }

        // TODO: Remove
        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }

            if(!(o instanceof RefMask))
            {
                return false;
            }

            final RefMask refMask = (RefMask) o;
            return PosStart == refMask.PosStart && PosEnd == refMask.PosEnd && BaseMask == refMask.BaseMask;
        }

        @Override
        public int hashCode()
        {
            int hash = PosStart;
            hash = PosEnd + 31 * hash;
            hash = BaseMask + 31 * hash;
            return hash;
        }
    }

    @VisibleForTesting
    public static class MergedHomopolymers
    {
        public final List<Homopolymer> RefHomopolymers;
        public final List<Homopolymer> ReadHomopolymers;

        private boolean mVariantInMergedHomopolymers;
        private List<RefMask> mRefMasks;

        public MergedHomopolymers(final List<Homopolymer> refHomopolymers, final List<Homopolymer> readHomopolymers,
                boolean variantInMergedHomopolymers)
        {
            RefHomopolymers = refHomopolymers;
            ReadHomopolymers = readHomopolymers;
            mVariantInMergedHomopolymers = variantInMergedHomopolymers;
            mRefMasks = Lists.newArrayList();
        }

        public MergedHomopolymers()
        {
            this(Lists.newArrayList(), Lists.newArrayList(), false);
        }

        public void setVariantInMergedHomopolymers()
        {
            mVariantInMergedHomopolymers = true;
        }

        public boolean variantInMergedHomopolymers()
        {
            return mVariantInMergedHomopolymers;
        }

        public void addRefMask(final RefMask refMask)
        {
            mRefMasks.add(refMask);
        }

        public List<RefMask> refMasks()
        {
            return mRefMasks;
        }

        // TODO: remove
        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof MergedHomopolymers))
            {
                return false;
            }
            final MergedHomopolymers that = (MergedHomopolymers) o;
            return mVariantInMergedHomopolymers == that.mVariantInMergedHomopolymers
                    && Objects.equals(RefHomopolymers, that.RefHomopolymers)
                    && Objects.equals(ReadHomopolymers, that.ReadHomopolymers) && Objects.equals(mRefMasks, that.mRefMasks);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(RefHomopolymers, ReadHomopolymers, mVariantInMergedHomopolymers, mRefMasks);
        }
    }

    @VisibleForTesting
    public static boolean isCleanSnv(final VariantReadContext readContext)
    {
        if(!readContext.variant().isSNV())
        {
            return false;
        }

        if(readContext.coreLength() != readContext.refBasesBytes().length)
        {
            return false;
        }

        for(int i = readContext.CoreIndexStart; i <= readContext.CoreIndexEnd; ++i)
        {
            if(i == readContext.VarIndex)
            {
                continue;  // the SNV base
            }

            byte readBase = readContext.readBasesBytes()[i];
            byte refBase = readContext.refBasesBytes()[i - readContext.CoreIndexStart];
            if(readBase != refBase)
            {
                return false;
            }
        }

        return true;
    }

    public static List<Integer> coreHomopolymerLengths(final VariantReadContext readContext)
    {
        List<Homopolymer> homopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        return homopolymers.stream().map(x -> x.Length).collect(Collectors.toList());
    }

}
