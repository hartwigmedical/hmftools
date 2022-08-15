package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;

import java.util.List;

import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.purple.ChromosomeArm;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.cn.SvCNData;

public class SvBreakend {

    private final SvVarData mSV;
    private int mChrPosIndex; // position in chromosome's breakend list
    private int mClusterChrPosIndex;
    private final String mChromosome;
    private final ChromosomeArm mArm ;
    private final String mChrArm;
    private final int mPosition;
    private final byte mOrientation;
    private boolean mUsesStart;
    private boolean mSglMapping;

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
        mSglMapping = false;
    }

    public SvBreakend(final SvVarData var, final SglMapping mapping)
    {
        mSV = var;
        mChrPosIndex = -1;
        mClusterChrPosIndex = -1;
        mUsesStart = false;
        mChromosome = mapping.Chromosome;
        mArm = getChromosomalArm(mapping.Chromosome, mapping.Position);
        mChrArm = makeChrArmStr(mChromosome, mArm);
        mPosition = mapping.Position;
        mOrientation = mapping.Orientation;
        mSglMapping = true;
    }

    public final SvVarData getSV() { return mSV; }
    public final SvCluster getCluster() { return mSV.getCluster(); }

    public final String chromosome() { return mChromosome; }
    public final ChromosomeArm arm() { return mArm; }
    public final String getChrArm() { return mChrArm; }

    public final int position() { return mPosition; }
    public final byte orientation() { return mOrientation; }
    public boolean usesStart() { return mUsesStart; }
    public StructuralVariantType type() { return mSV.type(); }
    public boolean isSglMapping() { return mSglMapping; }

    public final SvBreakend getOtherBreakend() { return mSV.getBreakend(!mUsesStart); }

    public void setChrPosIndex(int index) { mChrPosIndex = index; }
    public int getChrPosIndex() { return mChrPosIndex; }

    public void setClusterChrPosIndex(int index) { mClusterChrPosIndex = index; }
    public int getClusterChrPosIndex() { return mClusterChrPosIndex; }

    public static final String DIRECTION_CENTROMERE = "C";
    public static final String DIRECTION_TELOMERE = "T";

    public String direction() { return (mOrientation == 1) == (mArm == P_ARM) ? DIRECTION_TELOMERE : DIRECTION_CENTROMERE; }

    public final String toString()
    {
        return String.format("%s: %s %s:%d:%d", mSV.id(), mUsesStart ? "start" :"end", mChromosome, mPosition, mOrientation);
    }

    public int anchorDistance()
    {
        if(!mUsesStart && mSV.isSglBreakend())
            return 0;

        return mUsesStart ? mSV.getSvData().startAnchoringSupportDistance() : mSV.getSvData().endAnchoringSupportDistance();
    }

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
    public double jcn() { return mSV.jcn(); }
    public double jcnUncertainty() { return mSV.jcnUncertainty(); }

    public String homology()
    {
        if(!mUsesStart && mSV.isSglBreakend())
            return "";

        return mUsesStart ? mSV.getSvData().startHomologySequence() : mSV.getSvData().endHomologySequence();
    }

    public double getCopyNumber(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.CopyNumber;
    }

    public double majorAlleleJcn(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.majorAlleleJcn();
    }

    public double minorAlleleJcn(boolean usePrevious)
    {
        SvCNData cnData = mSV.getCopyNumberData(mUsesStart, usePrevious);

        if(cnData == null)
            return 0;

        return cnData.minorAlleleJcn();
    }

    public boolean isAssembledLink()
    {
        return mSV.hasAssemblyLink(mUsesStart);
    }

    public final DbPair getDBLink() { return mSV.getDBLink(mUsesStart); }
    public final List<LinkedPair> getLinkedPairs() { return mSV.getLinkedPairs(mUsesStart); }

    public boolean isFoldback()
    {
        return mSV.getFoldbackBreakend(mUsesStart) != null;
    }
    public final SvBreakend getFoldbackBreakend() { return mSV.getFoldbackBreakend(mUsesStart); }

    public boolean inLineElement()
    {
        return mSV.isLineElement(mUsesStart);
    }
    public boolean hasLineElement(LineElementType type) { return mSV.hasLineElement(type, mUsesStart); }

    public List<BreakendGeneData> getGenesList() { return mSV.getGenesList(mUsesStart); }

}
