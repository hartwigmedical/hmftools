package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.compareChromosomes;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class BasePosition implements Comparable<BasePosition>
{
    public final String Chromosome;
    public final int Position;

    public BasePosition(final String chromosome, final int position)
    {
        Chromosome = chromosome;
        Position = position;
    }

    @Override
    public int compareTo(final BasePosition other)
    {
        if(Chromosome.equals(other.Chromosome))
        {
            if(Position != other.Position)
                return Position < other.Position ? -1 : 1;

            return 0;
        }

        return compareChromosomes(Chromosome, other.Chromosome);
    }

    public boolean matches(final BasePosition other)
    {
        return Position == other.Position && Chromosome.equals(other.Chromosome);
    }
    public boolean matches(final String chromosome, final int position) { return Chromosome.equals(chromosome) && Position == position; }

    @Override
    public boolean equals(Object obj)
    {
        if(obj == this)
            return true;

        if(obj == null)
            return false;

        if(!getClass().equals(obj.getClass()))
            return false;

        BasePosition other = (BasePosition) obj;
        return matches(other);
    }

    @Override
    public int hashCode()
    {
        int result = 31 + Chromosome.hashCode();
        result = 31 * result + Position;
        return result;
    }

    public String toString() { return String.format("%s:%d", Chromosome, Position); }

    public static BasePosition from(final GenomePosition genomePosition) { return new BasePosition(genomePosition.chromosome(), genomePosition.position()); }
}
