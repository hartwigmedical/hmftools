package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;

import com.hartwig.hmftools.linx.cn.SvCNData;

public class SvBreakend {

    private final SvVarData mSV;
    private int mChrPosIndex; // position in chromosome's breakend list
    private int mClusterChrPosIndex;
    private final String mChromosome;
    private final String mArm ;
    private final String mChrArm;
    private final long mPosition;
    private final byte mOrientation;
    private boolean mUsesStart;

    public SvBreakend(final SvVarData var, boolean useStart)
    {
        mSV = var;
        mChrPosIndex = -1;
        mClusterChrPosIndex = -1;
        mUsesStart = useStart;
        mChromosome = var.chromosome(useStart);
        mArm = var.arm(useStart);
        mChrArm = makeChrArmStr(mChromosome, mArm);
        mPosition = var.position(useStart);
        mOrientation = var.orientation(useStart);
    }

    public final SvVarData getSV() { return mSV; }
    public final SvCluster getCluster() { return mSV.getCluster(); }

    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public final String getChrArm() { return mChrArm; }

    public final long position() { return mPosition; }
    public final byte orientation() { return mOrientation; }
    public boolean usesStart() { return mUsesStart; }

    public final SvBreakend getOtherBreakend() { return mSV.getBreakend(!mUsesStart); }

    public void setChrPosIndex(int index) { mChrPosIndex = index; }
    public int getChrPosIndex() { return mChrPosIndex; }

    public void setClusterChrPosIndex(int index) { mClusterChrPosIndex = index; }
    public int getClusterChrPosIndex() { return mClusterChrPosIndex; }

    public static String DIRECTION_CENTROMERE = "C";
    public static String DIRECTION_TELOMERE = "T";

    public String direction() { return (mOrientation == 1) == (mArm == CHROMOSOME_ARM_P) ? DIRECTION_TELOMERE : DIRECTION_CENTROMERE; }

    public final String toString()
    {
        return mSV.posId(mUsesStart);
    }

    public int getMinTemplatedLength() { return mSV.getMinTemplatedLength(mUsesStart); }

    public double copyNumberLowSide()
    {
        return mSV.copyNumber(mUsesStart) - mSV.copyNumberChange(mUsesStart);
    }

    // for convenience
    public double copyNumberChange()
    {
        return mSV.copyNumberChange(mUsesStart);
    }
    public double copyNumber()
    {
        return mSV.copyNumber(mUsesStart);
    }
    public double ploidy() { return mSV.ploidy(); }
    public double ploidyUncertainty() { return mSV.ploidyUncertainty(); }

    public double getCopyNumber(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.CopyNumber;
    }

    public double majorAllelePloidy(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.majorAllelePloidy();
    }

    public double minorAllelePloidy(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.minorAllelePloidy();
    }

    public boolean isAssembledLink()
    {
        return mSV.hasAssemblyLink(mUsesStart);
    }

    public final SvLinkedPair getDBLink() { return mSV.getDBLink(mUsesStart); }

    public boolean isFoldback()
    {
        return mSV.getFoldbackBreakend(mUsesStart) != null;
    }
    public final SvBreakend getFoldbackBreakend() { return mSV.getFoldbackBreakend(mUsesStart); }

}
