package com.hartwig.hmftools.common.utils.pcf;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class PCFPosition implements GenomePosition
{
    public final PCFSource Source;
    public final String Chromosome;
    public final int Position;

    private int mMinPosition;
    private int mMaxPosition;

    public PCFPosition(final PCFSource source, final String chromosome, final int position)
    {
        this(source, chromosome, position, position, position);
    }

    public PCFPosition(final PCFSource source, final String chromosome, final int position, final int minPosition, final int maxPosition)
    {
        Source = source;
        Chromosome = chromosome;
        Position = position;
        mMinPosition = minPosition;
        mMaxPosition = maxPosition;
    }

    public int minPosition() { return mMinPosition; }
    public int maxPosition() { return mMaxPosition; }

    public void setMinPosition(int position) { mMinPosition = position; }
    public void setMaxPosition(int position) { mMaxPosition = position; }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public int position()
    {
        return Position;
    }

    public boolean isSegmentEnd() { return Position == minPosition(); }
    public boolean isSegmentStart() { return !isSegmentEnd(); }

    public String toString()
    {
        return String.format("src(%s) pos(%s:%d) range(%d - %d) segment %s",
                Source, Chromosome, Position, mMinPosition, mMaxPosition, Position == mMinPosition ? "end" : "start");
    }
}
