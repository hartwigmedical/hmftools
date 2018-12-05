package com.hartwig.hmftools.svanalysis.types;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;

public class SvBreakend {

    private final SvVarData mSV;
    private int mChrPosIndex; // position in chromosome
    private final String mChromosome;
    private final String mArm ;
    private final String mChrArm;
    private final Long mPosition;
    private final Byte mOrientation;
    private boolean mUsesStart;

    public SvBreakend(final SvVarData var, boolean useStart)
    {
        mSV = var;
        mChrPosIndex = -1;
        mUsesStart = useStart;
        mChromosome = var.chromosome(useStart);
        mArm = var.arm(useStart);
        mChrArm = makeChrArmStr(mChromosome, mArm);
        mPosition = var.position(useStart);
        mOrientation = var.orientation(useStart);
    }

    public final SvVarData getSV() { return mSV; }
    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public final String getChrArm() { return mChrArm; }

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
