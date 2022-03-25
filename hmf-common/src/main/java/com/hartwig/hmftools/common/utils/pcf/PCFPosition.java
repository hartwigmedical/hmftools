package com.hartwig.hmftools.common.utils.pcf;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class PCFPosition implements GenomePosition {

    public final PCFSource Source;
    public final String Chromosome;
    public final int Position;

    private int mMinPosition;
    private int mMaxPosition;

    public PCFPosition(final PCFSource source, final String chromosome, final int position)
    {
        Source = source;
        Chromosome = chromosome;
        Position = position;
        mMinPosition = position;
        mMaxPosition = position;
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

    public String toString()
    {
        return String.format("src(%s) pos(%s:%d) range(%d - %d)", Source, Chromosome, Position, mMinPosition, mMaxPosition);
    }
}
