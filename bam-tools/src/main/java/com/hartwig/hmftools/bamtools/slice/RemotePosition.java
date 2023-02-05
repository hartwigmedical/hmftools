package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import htsjdk.samtools.SAMRecord;

public class RemotePosition
{
    public final String ReadId;
    public final String Chromosome;
    public final int Position;

    public RemotePosition(final String readId, final String chromosome, final int position)
    {
        ReadId = readId;
        Chromosome = chromosome;
        Position = position;
    }

    public boolean positionsMatch(final RemotePosition other)
    {
        return Chromosome.equals(other.Chromosome) && Position == other.Position;
    }

    public boolean positionsMatch(final SAMRecord read)
    {
        return Chromosome.equals(read.getReferenceName()) && Position == read.getAlignmentStart();
    }

    public String toString()
    {
        return format("%s:%d %s", Chromosome, Position, ReadId);
    }

}
