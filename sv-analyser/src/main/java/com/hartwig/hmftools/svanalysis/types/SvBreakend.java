package com.hartwig.hmftools.svanalysis.types;

import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.MIN_LOH_CN;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;

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
    public final SvVarData getOrigSV() { return mSV.getOrigSV(); }
    public final SvBreakend getOrigBreakend() { return getOrigSV().getBreakend(mUsesStart); }
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

    public double getCopyNumber(boolean lowSide)
    {
        if((mOrientation == 1 && lowSide) || (mOrientation == -1 && !lowSide))
        {
            return mSV.copyNumber(mUsesStart);
        }
        else
        {
            return mSV.copyNumber(mUsesStart) - mSV.copyNumberChange(mUsesStart);
        }
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

    public double actualBaf(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.ActualBaf;
    }

    public boolean isAssembledLink()
    {
        return mSV.getAssemblyMatchType(mUsesStart) == ASSEMBLY_MATCH_MATCHED;
    }


}
