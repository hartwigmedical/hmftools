package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.CigarUtils.cigarFromStr;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;

import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.samtools.ClippedSide;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.sv.Direction;

import htsjdk.samtools.Cigar;

public class RemoteJunction
{
    public final String Chromosome;
    public final int Position;
    public final Direction Orientation;
    public int Fragments;

    public RemoteJunction(final String chromosome, final int position, final Direction orientation)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        Fragments = 1;
    }

    public static RemoteJunction fromSupplementaryData(final SupplementaryReadData suppData)
    {
        // use the cigar if the alignment suggest the soft-clip is on the right side
        Cigar remoteCigar = cigarFromStr(suppData.Cigar);
        ClippedSide scSide = ClippedSide.fromCigar(remoteCigar, true);

        if(scSide.isLeft())
        {
            return new RemoteJunction(suppData.Chromosome, suppData.Position, Direction.REVERSE);
        }
        else
        {
            int skippedBases = remoteCigar.getCigarElements().stream()
                    .filter(x -> x.getOperator() == N || x.getOperator() == M || x.getOperator() == D)
                    .mapToInt(x -> x.getLength()).sum();

            int remoteJunctionPos = suppData.Position + skippedBases - 1;
            return new RemoteJunction(suppData.Chromosome, remoteJunctionPos, Direction.FORWARDS);
        }
    }

    public boolean matches(final RemoteJunction other)
    {
        return Chromosome.equals(other.Chromosome) && Position == other.Position && Orientation == other.Orientation;
    }

    public boolean matches(final String chromosome, final int position, final Direction orientation)
    {
        return Chromosome.equals(chromosome) && Position == position && Orientation == orientation;
    }

    public static void addRemoteJunction(final List<RemoteJunction> remoteJunctions, final RemoteJunction remoteJunction)
    {
        if(remoteJunctions.stream().noneMatch(x -> x.matches(remoteJunction)))
            remoteJunctions.add(remoteJunction);
    }

    public String toString() { return format("loc(%s:%d:%d) reads(%d)", Chromosome, Position, Orientation.Step, Fragments); }

    public static class RemoteJunctionSorter implements Comparator<RemoteJunction>
    {
        // sorts by support descending
        public int compare(final RemoteJunction first, final RemoteJunction second)
        {
            if(first.Fragments != second.Fragments)
                return first.Fragments < second.Fragments ? 1 : -1;

            return 0;
        }
    }
}
