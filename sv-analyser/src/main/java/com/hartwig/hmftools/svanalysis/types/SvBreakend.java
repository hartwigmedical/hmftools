package com.hartwig.hmftools.svanalysis.types;

public class SvBreakend {

    private final SvVarData mSV;
    private final String mChromosome;
    private final Long mPosition;
    private final Byte mOrientation;
    private boolean mUsesStart;
    private int mCount;

    public SvBreakend(final SvVarData var, boolean useStart)
    {
        mSV = var;
        mUsesStart = useStart;
        mChromosome = var.chromosome(useStart);
        mPosition = var.position(useStart);
        mOrientation = var.orientation(useStart);
        mCount = 1;
    }

    public SvBreakend(final String chromosome, final Long position, final Byte orientation)
    {
        mSV = null;
        mUsesStart = false;
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
        mCount = 1;
    }

    public final SvVarData getSV() { return mSV; }
    public String chromosome() { return mChromosome; }
    public final Long position() { return mPosition; }
    public final Byte orientation() { return mOrientation; }
    public boolean usesStart() { return mUsesStart; }

    public int getCount() { return mCount; }
    public void addToCount(int change) { mCount += change; }

    public final String toString()
    {
        if(mSV == null)
            return String.format("%s %d:%d", mChromosome, mOrientation, mPosition);

        return mSV.posId(mUsesStart);
    }
}
