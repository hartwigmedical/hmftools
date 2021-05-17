package com.hartwig.hmftools.lilac.read;

import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class Indel
{
    public final boolean IsInsert;
    public final boolean IsDelete;
    public final int Length;
    public final String Contig;
    public final int Position;
    public final String Ref;
    public final String Alt;

    private final int mHashcode;

    public Indel(final String contig, int position, final String ref, final String alt)
    {
        Contig = contig;
        Position = position;
        Ref = ref;
        Alt = alt;
        IsInsert = Ref.length() < Alt.length();
        IsDelete = !IsInsert;
        Length = Alt.length() - Ref.length();

        mHashcode = toString().hashCode();
    }

    public static Indel fromString(final String line)
    {
        final String[] items = line.split("\t");
        return new Indel(items[0], Integer.parseInt(items[1]), items[2], items[3]);
    }

    public String toString()
    {
        return Contig + ':' + Position + ' ' + Ref + '>' + Alt;
    }

    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if (!(other instanceof Indel))
            return false;

        return hashCode() == other.hashCode();
    }

    public int hashCode() { return mHashcode; }

    public final boolean match(final Indel other)
    {
        if(equals(other))
            return true;

        if(!Contig.equals(other.Contig))
            return false;

        if(Position != other.Position)
            return false;

        if(Ref.length() != other.Ref.length())
            return false;

        if(Alt.length() != other.Alt.length())
            return false;

        if(Alt.length() > Ref.length())
        {
            // check the insert
            if(!Alt.substring(1).equals(other.Alt.substring(1)))
                return false;
        }

        return true;
    }
}
