package com.hartwig.hmftools.svassembly.models;

public class Alignment
{
    public final String Chromosome;
    public final int ReferenceStartPosition;
    public final int SequenceStartPosition;
    public final int Length;
    public final boolean Inverted;
    public final int Quality;

    public Alignment(final String chromosome, final int referenceStartPosition,
            final int sequenceStartPosition, final int length, final boolean inverted,
            final int quality)
    {
        assert sequenceStartPosition > 0;
        Chromosome = chromosome;
        ReferenceStartPosition = referenceStartPosition;
        SequenceStartPosition = sequenceStartPosition;
        Length = length;
        Inverted = inverted;
        Quality = quality;
    }

    public static Alignment unmapped(final int length)
    {
        return unmapped(1, length);
    }

    public static Alignment unmapped(final int startPosition, final int length)
    {
        return new Alignment("?", 0, startPosition, length, false, 0);
    }

    public static Alignment insert(final int startPosition, final int length)
    {
        return new Alignment("*", 0, startPosition, length, false, 0);
    }

    public static Alignment gap(final int startPosition)
    {
        return new Alignment("-", 0, startPosition, 1, false, 0);
    }

    public boolean isUnmapped()
    {
        return Chromosome.equals("*") || Chromosome.equals("?") || Chromosome.equals("-");
    }

    public boolean isMapped()
    {
        return !isUnmapped();
    }

    public boolean includes(final String chromosome, final int position)
    {
        return Chromosome.equals(chromosome) && ReferenceStartPosition <= position && position < ReferenceStartPosition + Length;
    }

    public boolean includes(final String chromosome, final int position, final int length)
    {
        if(!Chromosome.equals(chromosome))
            return false;

        return ReferenceStartPosition < (position + length) && (ReferenceStartPosition + Length) > position;
    }

    @Override
    public String toString()
    {
        return String.format("%s:%d (%d bp%s) @ %d", Chromosome, ReferenceStartPosition, Length, Inverted
                ? ", inverted"
                : "", SequenceStartPosition);
    }
}
