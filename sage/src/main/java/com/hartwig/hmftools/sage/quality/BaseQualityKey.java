package com.hartwig.hmftools.sage.quality;

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

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof BaseQualityKey))
            return false;

        BaseQualityKey otherKey = (BaseQualityKey)other;
        return matches(otherKey.Ref, otherKey.Alt, otherKey.Quality, otherKey.TrinucleotideContext);
        // return hashCode() == other.hashCode();
    }

    public boolean matches(byte ref, byte alt, byte quality, byte[] trinucleotideContext)
    {
        if(Ref != ref || Alt != alt || Quality != quality)
            return false;

        return (trinucleotideContext[0] == TrinucleotideContext[0] && trinucleotideContext[1] == TrinucleotideContext[1]
                &&  trinucleotideContext[2] == TrinucleotideContext[2]);
    }

    public String toString() { return String.format("var(%c->%c) cxt(%s) qual(%d)",
            (char)Ref, (char)Alt, TrinucleotideContext != null ? new String(TrinucleotideContext) : "", (int)Quality);}
}
