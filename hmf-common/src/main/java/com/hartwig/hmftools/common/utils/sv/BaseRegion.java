package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.NoSuchElementException;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class BaseRegion implements Cloneable, Comparable<BaseRegion>
{
    private int mStart;
    private int mEnd;

    public BaseRegion(final int posStart, final int posEnd)
    {
        mStart = posStart;
        mEnd = posEnd;
    }

    public static BaseRegion from(final GenomeRegion region) { return new BaseRegion(region.start(), region.end()); }

    public int start() { return mStart; }
    public int end() { return mEnd; }

    public int position(int which)
    {
        if(which == SE_START)
        {
            return mStart;
        }
        else if(which == SE_END)
        {
            return mEnd;
        }
        throw new NoSuchElementException();
    }

    public void setStart(int pos) { mStart = pos; }
    public void setEnd(int pos) { mEnd = pos; }

    public int baseLength() { return length() + 1; }
    public int length() { return mEnd - mStart; }

    public boolean hasValidPositions() { return mStart > 0 & mEnd >= mStart; }

    public boolean overlaps(final BaseRegion other)
    {
        return positionsOverlap(mStart, mEnd, other.mStart, other.mEnd);
    }

    public boolean overlaps(final ChrBaseRegion other)
    {
        // assumes chromosome check is not relevant
        return positionsOverlap(mStart, mEnd, other.start(), other.end());
    }

    public boolean containsPosition(int position) { return positionWithin(position, start(), end()); }

    public boolean matches(final BaseRegion other)
    {
        return start() == other.start() && end() == other.end();
    }

    public String toString() { return String.format("%d-%d", mStart, mEnd); }

    @Override
    public BaseRegion clone()
    {
        try
        {
            BaseRegion br = (BaseRegion) super.clone();
            br.mStart = mStart;
            br.mEnd = mEnd;
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
    public int hashCode()
    {
        int result = 31 + mStart;
        result = 31 * result + mEnd;
        return result;
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

