package com.hartwig.hmftools.bamtools.slice;

import static java.lang.String.format;

import com.hartwig.hmftools.common.sage.GeneDepth;

import htsjdk.samtools.SAMRecord;

public class RemotePosition implements Comparable<RemotePosition>
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

    @Override
    public int compareTo(final RemotePosition other)
    {
        if(Position == other.Position)
            return 0;

        return Position < other.Position ? -1 : 1;
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
