package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class BaseRegion implements Cloneable, Comparable<BaseRegion>
{
    public int[] Positions;

    public BaseRegion(final int[] positions)
    {
        Positions = positions;
    }

    public BaseRegion(final int posStart, final int posEnd)
    {
        Positions = new int[] { posStart, posEnd };
    }

    public static BaseRegion from(final GenomeRegion region) { return new BaseRegion((int)region.start(), (int)region.end()); }

    public int start() { return Positions[SE_START]; }
    public int end() { return Positions[SE_END]; }

    public void setPosition(int position, int index) { Positions[index] = position; }
    public void setStart(int pos) { setPosition(pos, SE_START); }
    public void setEnd(int pos) { setPosition(pos, SE_END); }

    public int baseLength() { return length() + 1; }
    public int length() { return Positions[SE_END] - Positions[SE_START]; }

    public boolean hasValidPositions() { return Positions[SE_START] > 0 & Positions[SE_END] >= Positions[SE_START]; }

    public boolean overlaps(final BaseRegion other)
    {
        return positionsOverlap(Positions[SE_START], Positions[SE_END], other.Positions[SE_START], other.Positions[SE_END]);
    }

    public boolean containsPosition(int position) { return positionWithin(position, start(), end()); }

    public boolean matches(final BaseRegion other)
    {
        return start() == other.start() && end() == other.end();
    }

    public String toString() { return String.format("%d-%d", Positions[SE_START], Positions[SE_END]); }

    @Override
    public Object clone()
    {
        try
        {
            BaseRegion br = (BaseRegion) super.clone();
            br.Positions = Positions.clone();
            return br;
        }
        catch (CloneNotSupportedException e)
        {
            // Will not happen in this case
            return null;
        }
    }

    @Override
    public boolean equals(Object obj)
    {
        // same instance
        if (obj == this) { return true; }
        // null
        if (obj == null) { return false; }
        // type
        if (!getClass().equals(obj.getClass())) { return false; }
        // cast and compare state
        BaseRegion other = (BaseRegion) obj;
        return matches(other);
    }

    @Override
    public int compareTo(@NotNull final BaseRegion other)
    {
        if (start() < other.start())
        {
            return -1;
        }
        else if (start() == other.start())
        {
            return 0;
        }
        return 1;
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
        return (innerStart <= innerEnd && innerStart >= outerStart && innerEnd <= outerEnd);
    }

}

