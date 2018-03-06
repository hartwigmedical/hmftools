package com.hartwig.hmftools.svannotation.analysis;

public class SvBreakend {

    private final String mChromosome;
    private final Long mPosition;
    private final Byte mOrientation;
    private int mCount;

    public SvBreakend(final String chromosome, final Long position, final Byte orientation)
    {
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
        mCount = 1;
    }

    public String chromosome() { return mChromosome; }
    public final Long position() { return mPosition; }
    public final Byte orientation() { return mOrientation; }

    public int getCount() { return mCount; }
    public void addToCount(int change) { mCount += change; }
}
