package com.hartwig.hmftools.sage.quality;

import java.util.Arrays;

import com.google.common.primitives.Bytes;

public class QualityRecalibrationKey
{
    public final byte Ref;
    public final byte Alt;
    public final byte[] TrinucleotideContext;
    public final byte Quality;

    private int mHashCode;

    public QualityRecalibrationKey(final byte ref, final byte alt, final byte[] trinucleotideContext, final byte quality)
    {
        Ref = ref;
        Alt = alt;
        TrinucleotideContext = trinucleotideContext;
        Quality = quality;

        mHashCode = calcHashCode();
    }

    public static QualityRecalibrationKey from(final QualityCounterKey key)
    {
        return new QualityRecalibrationKey(key.Ref, key.Alt, key.TrinucleotideContext, key.Quality);
    }

    public int hashCode() { return mHashCode; }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof QualityRecalibrationKey))
            return false;

        return hashCode() == other.hashCode();
    }

    public String toString() { return String.format("var(%c->%c) cxt(%s) qual(%d)",
            (char)Ref, (char)Alt, new String(TrinucleotideContext), Quality);}

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
