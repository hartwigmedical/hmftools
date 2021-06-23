package com.hartwig.hmftools.sage.quality;

import java.util.Arrays;

import com.google.common.primitives.Bytes;

public class QualityCounterKey implements Comparable<QualityCounterKey>
{
    public final byte Ref;
    public final byte Alt;
    public final byte Quality;
    public final int Position;

    public final byte[] TrinucleotideContext;

    private int mHashCode;

    public QualityCounterKey(final byte ref, final byte alt, final byte qual, final int position, final byte[] trinucleotideContext)
    {
        Ref = ref;
        Alt = alt;
        Quality = qual;
        Position = position;
        TrinucleotideContext = trinucleotideContext;

        mHashCode = calcHashCode();
    }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof QualityCounterKey))
            return false;

        return hashCode() == other.hashCode();
    }

    @Override
    public int compareTo(final QualityCounterKey other)
    {
        return compareTo(other, true);
    }

    public int compareTo(final QualityCounterKey other, boolean compareAll)
    {
        if(hashCode() == other.hashCode())
            return 0;

        int refCompare = Byte.compare(Ref, other.Ref);
        if(refCompare != 0)
            return refCompare;

        int altCompare = Byte.compare(Alt, other.Alt);
        if(altCompare != 0)
            return altCompare;

        if(compareAll)
        {
            int pos = Integer.compare(Position, other.Position);
            if(pos != 0)
                return pos;

            int qualCompare = Byte.compare(Quality, other.Quality);
            if(qualCompare != 0)
                return qualCompare;
        }

        if(TrinucleotideContext == null || TrinucleotideContext.length < 3
        || other.TrinucleotideContext == null || other.TrinucleotideContext.length < 3)
        {
            return 0;
        }

        int triOne = Byte.compare(TrinucleotideContext[0], other.TrinucleotideContext[0]);
        if(triOne != 0)
            return triOne;

        int triTwo = Byte.compare(TrinucleotideContext[1], other.TrinucleotideContext[1]);
        if(triTwo != 0)
            return triTwo;

        int triThree = Byte.compare(TrinucleotideContext[2], other.TrinucleotideContext[2]);
        if(triThree != 0)
            return triThree;

        return 0;
    }

    public int hashCode() { return mHashCode; }

    public String toString() { return String.format("var(%c->%c) pos(%d) cxt(%s) qual(%d)",
            (char)Ref, (char)Alt, Position, new String(TrinucleotideContext), (int)Quality);}

    private int calcHashCode()
    {
        String hash = String.format("%c%c%d%s%d",
                (char)Ref, (char)Alt, Position, TrinucleotideContext != null ? new String(TrinucleotideContext) : "", (int)Quality);

        return hash.hashCode();
        /*
        int h = 5381;
        h += (h << 5) + Bytes.hashCode(Ref);
        h += (h << 5) + Bytes.hashCode(Alt);
        h += (h << 5) + Bytes.hashCode(Quality);
        h += (h << 5) + Position;
        h += (h << 5) + Arrays.hashCode(TrinucleotideContext);
        return h;
        */
    }
}
