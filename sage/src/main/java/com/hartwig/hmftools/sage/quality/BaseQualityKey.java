package com.hartwig.hmftools.sage.quality;

import java.util.Arrays;
import java.util.Objects;

public class BaseQualityKey
{
    public final byte Ref;
    public final byte Alt;
    public final byte[] TrinucleotideContext;
    public final byte Quality;

    public BaseQualityKey(final byte ref, final byte alt, final byte[] trinucleotideContext, final byte quality)
    {
        Ref = ref;
        Alt = alt;
        TrinucleotideContext = trinucleotideContext;
        Quality = quality;
    }

    @Override
    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof BaseQualityKey))
            return false;

        BaseQualityKey otherKey = (BaseQualityKey)other;
        return matches(otherKey.Ref, otherKey.Alt, otherKey.Quality, otherKey.TrinucleotideContext);
    }

    public boolean matches(byte ref, byte alt, byte quality, byte[] trinucleotideContext)
    {
        if(Ref != ref || Alt != alt || Quality != quality)
            return false;

        return Arrays.equals(trinucleotideContext, TrinucleotideContext);
    }

    @Override
    public int hashCode()
    {
        int result = Objects.hash(Ref, Alt, Quality);
        result = 31 * result + Arrays.hashCode(TrinucleotideContext);
        return result;
    }

    public String toString() { return String.format("var(%c->%c) cxt(%s) qual(%d)",
            (char)Ref, (char)Alt, TrinucleotideContext != null ? new String(TrinucleotideContext) : "", (int)Quality);}
}
