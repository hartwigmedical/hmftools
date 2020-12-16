package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class SvRegion
{
    public final String Chromosome;
    public final int[] Positions;

    public SvRegion(final String chromosome, final int[] positions)
    {
        Chromosome = chromosome;
        Positions = positions;
    }

    public SvRegion(final String chromosome, final int posStart, final int posEnd)
    {
        Chromosome = chromosome;
        Positions = new int[] { posStart, posEnd };
    }

    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }

    public void setPosition(int position, int index) { Positions[index] = position; }
    public void setStart(int pos) { setPosition(pos, SE_START); }
    public void setEnd(int pos) { setPosition(pos, SE_END); }

    public int baseLength() { return length() + 1; }
    public int length() { return Positions[SE_END] - Positions[SE_START]; }

    public boolean isValid() { return HumanChromosome.contains(Chromosome) && hasValidPositions(); }
    public boolean hasValidPositions() { return Positions[SE_START] > 0 & Positions[SE_END] >= Positions[SE_START]; }

    public boolean overlaps(final SvRegion other)
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

    public String toString() { return String.format("%s:%d-%d", Chromosome, Positions[SE_START], Positions[SE_END]); }

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
