package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class BaseRegion implements Comparable<BaseRegion>
{
    public final String Chromosome;
    public final int[] Positions;

    public BaseRegion(final String chromosome, final int[] positions)
    {
        Chromosome = chromosome;
        Positions = positions;
    }

    public BaseRegion(final String chromosome, final int posStart, final int posEnd)
    {
        Chromosome = chromosome;
        Positions = new int[] { posStart, posEnd };
    }

    public static BaseRegion from(final GenomeRegion region) { return new BaseRegion(region.chromosome(), (int)region.start(), (int)region.end()); }

    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }
    public String chromosome() { return Chromosome; }

    public void setPosition(int position, int index) { Positions[index] = position; }
    public void setStart(int pos) { setPosition(pos, SE_START); }
    public void setEnd(int pos) { setPosition(pos, SE_END); }

    public int baseLength() { return length() + 1; }
    public int length() { return Positions[SE_END] - Positions[SE_START]; }

    public boolean isValid() { return HumanChromosome.contains(Chromosome) && hasValidPositions(); }
    public boolean hasValidPositions() { return Positions[SE_START] > 0 & Positions[SE_END] >= Positions[SE_START]; }

    public boolean overlaps(final BaseRegion other)
    {
        if(!Chromosome.equals(other.Chromosome))
            return false;

        return positionsOverlap(Positions[SE_START], Positions[SE_END], other.Positions[SE_START], other.Positions[SE_END]);
    }

    public boolean containsPosition(int position) { return positionWithin(position, start(), end()); }

    public boolean containsPosition(final String chromosome, int position)
    {
        return Chromosome.equals(chromosome) && positionWithin(position, start(), end());
    }

    public boolean matches(final BaseRegion other)
    {
        return Chromosome.equals(other.Chromosome) && start() == other.start() && end() == other.end();
    }

    public String toString() { return String.format("%s:%d-%d", Chromosome, Positions[SE_START], Positions[SE_END]); }

    @Override
    public int compareTo(@NotNull final BaseRegion other)
    {
        if(Chromosome.equals(other.Chromosome))
        {
            if (start() < other.start())
            {
                return -1;
            } else if (start() == other.start()) {
                return 0;
            }
            return 1;
        }

        return ContigComparator.INSTANCE.compare(Chromosome, other.Chromosome);
    }

    // utility methods relating to position comparisons
    public static boolean positionsOverlap(int posStart1, int posEnd1, int posStart2, int posEnd2)
    {
        return !(posStart1 > posEnd2 || posEnd1 < posStart2);
    }

    public static boolean positionWithin(int position, int otherPosStart, int otherPosEnd)
    {
        return (position >= otherPosStart && position <= otherPosEnd);
    }

    public static boolean positionsWithin(int innerStart, int innerEnd, int outerStart, int outerEnd)
    {
        return (innerStart < innerEnd && innerStart >= outerStart && innerEnd <= outerEnd);
    }

}
