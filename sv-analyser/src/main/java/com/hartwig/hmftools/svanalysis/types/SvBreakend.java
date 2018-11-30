package com.hartwig.hmftools.svanalysis.types;

public class SvBreakend {

    private final SvVarData mSV;
    private int mChrPosIndex; // position in chromosome
    private final String mChromosome;
    private final Long mPosition;
    private final Byte mOrientation;
    private boolean mUsesStart;

    public SvBreakend(final SvVarData var, boolean useStart)
    {
        mSV = var;
        mChrPosIndex = -1;
        mUsesStart = useStart;
        mChromosome = var.chromosome(useStart);
        mPosition = var.position(useStart);
        mOrientation = var.orientation(useStart);
    }

    public final SvVarData getSV() { return mSV; }
    public String chromosome() { return mChromosome; }
    public final Long position() { return mPosition; }
    public final Byte orientation() { return mOrientation; }
    public boolean usesStart() { return mUsesStart; }

    public void setChrPosIndex(int index) { mChrPosIndex = index; }
    public int getChrPosIndex() { return mChrPosIndex; }

    public final String toString()
    {
        if(mSV == null)
            return String.format("%s %d:%d", mChromosome, mOrientation, mPosition);

        return mSV.posId(mUsesStart);
    }
}
