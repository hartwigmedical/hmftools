package com.hartwig.hmftools.lilac.read;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;

import java.util.Objects;

public class Indel
{
    private static final int INVALID_READ_INDEX = -1;

    public final boolean IsInsert;
    public final boolean IsDelete;
    public final int Length;
    public final String Contig;
    public final int Position;
    public final int ReadIndex;
    public final String Ref;
    public final String Alt;
    public final boolean Ignore;

    public Indel(final String contig, int position, int readIndex, final String ref, final String alt, boolean ignore)
    {
        Contig = contig;
        Position = position;
        ReadIndex = readIndex;
        Ref = ref;
        Alt = alt;
        IsInsert = Ref.length() < Alt.length();
        IsDelete = !IsInsert;
        Length = Alt.length() - Ref.length();
        Ignore = ignore;
    }

    public Indel(final String contig, int position, final String ref, final String alt)
    {
        this(contig, position, INVALID_READ_INDEX, ref, alt, false);
    }

    public static Indel fromString(final String line)
    {
        final String[] items = line.split(CSV_DELIM);
        return new Indel(items[0], Integer.parseInt(items[1]), items[2], items[3]);
    }

    @Override
    public String toString()
    {
        return Contig + ':' + Position + ' ' + Ref + '>' + Alt;
    }

    @Override
    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof Indel))
            return false;

        return toString().equals(other.toString()) && Ignore == ((Indel) other).Ignore;
    }

    @Override
    public int hashCode() { return Objects.hash(toString(), Ignore); }

    public final boolean match(final Indel other)
    {
        if(equals(other))
            return true;

        if(Ignore != other.Ignore)
            return false;

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
