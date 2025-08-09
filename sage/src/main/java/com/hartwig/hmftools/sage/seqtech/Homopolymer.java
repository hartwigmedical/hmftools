package com.hartwig.hmftools.sage.seqtech;

import java.util.List;

import com.google.common.collect.Lists;

public class Homopolymer
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
}
