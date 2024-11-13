package com.hartwig.hmftools.common.qual;

import ngs.Read;

import java.util.Arrays;
import java.util.Objects;

public class BqrKey
{
    public final byte Ref;
    public final byte Alt;
    public final byte[] TrinucleotideContext;
    public final byte Quality;
    public final BqrReadType ReadType;
    public final BqrReadStrand ReadStrand;

    public BqrKey(final byte ref, final byte alt, final byte[] trinucleotideContext, final byte quality, final BqrReadType readType, final BqrReadStrand readStrand)
    {
        Ref = ref;
        Alt = alt;
        TrinucleotideContext = trinucleotideContext;
        Quality = quality;
        ReadType = readType;
        ReadStrand = readStrand;
    }

    @Override
    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof BqrKey))
            return false;

        BqrKey otherKey = (BqrKey)other;
        return matches(otherKey.Ref, otherKey.Alt, otherKey.Quality, otherKey.TrinucleotideContext, otherKey.ReadType, otherKey.ReadStrand);
    }

    public boolean matches(byte ref, byte alt, byte quality, byte[] trinucleotideContext, BqrReadType readType, BqrReadStrand readStrand)
    {
        if(Ref != ref || Alt != alt || Quality != quality || ReadType != readType || ReadStrand != readStrand)
            return false;

        return Arrays.equals(trinucleotideContext, TrinucleotideContext);
    }

    @Override
    public int hashCode()
    {
        int result = Objects.hash(Ref, Alt, Quality, ReadType.ordinal(), ReadStrand.ordinal());
        result = 31 * result + Arrays.hashCode(TrinucleotideContext);
        return result;
    }

    public String toString()
    {
        return String.format("var(%c->%c) cxt(%s) qual(%d) type(%s) strand(%s)",
            (char)Ref, (char)Alt, TrinucleotideContext != null ? new String(TrinucleotideContext) : "", (int)Quality, ReadType, ReadStrand);
    }
}
