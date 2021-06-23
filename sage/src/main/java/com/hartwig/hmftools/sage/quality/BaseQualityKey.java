package com.hartwig.hmftools.sage.quality;

import java.util.Arrays;

import com.google.common.primitives.Bytes;

public class BaseQualityKey
{
    public final byte Ref;
    public final byte Alt;
    public final byte[] TrinucleotideContext;
    public final byte Quality;

    private int mHashCode;

    public BaseQualityKey(final byte ref, final byte alt, final byte[] trinucleotideContext, final byte quality)
    {
        Ref = ref;
        Alt = alt;
        TrinucleotideContext = trinucleotideContext;
        Quality = quality;

        mHashCode = calcHashCode();
    }

    public int hashCode() { return mHashCode; }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof BaseQualityKey))
            return false;

        return hashCode() == other.hashCode();
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

    private int calcHashCode()
    {
        int h = 5381;
        h += (h << 5) + Bytes.hashCode(Ref);
        h += (h << 5) + Bytes.hashCode(Alt);
        h += (h << 5) + Bytes.hashCode(Quality);
        h += (h << 5) + Arrays.hashCode(TrinucleotideContext);
        return h;
    }

}
