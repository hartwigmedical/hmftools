package com.hartwig.hmftools.purple.copynumber.sv;

import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

public class StructuralVariantLegCopyNumber implements StructuralVariantLeg
{
    private final StructuralVariantLeg mLeg;

    private Double mLeftCopyNumber;
    private Double mRightCopyNumber;

    public StructuralVariantLegCopyNumber(final StructuralVariantLeg leg, final Double leftCopyNumber, final Double rightCopyNumber)
    {
        mLeg = leg;
        mLeftCopyNumber = leftCopyNumber;
        mRightCopyNumber = rightCopyNumber;
    }

    public double leftCopyNumberOrZero() { return mLeftCopyNumber != null ? mLeftCopyNumber : 0; }
    public Double leftCopyNumber() { return mLeftCopyNumber; }
    public void setLeftCopyNumber(final Double value) { mLeftCopyNumber = value; }

    public double rightCopyNumberOrZero() { return mRightCopyNumber != null ? mRightCopyNumber : 0; }
    public Double rightCopyNumber() { return mRightCopyNumber; }
    public void setRightCopyNumber(final Double value) { mRightCopyNumber = value; }

    @Override
    public String chromosome() { return mLeg.chromosome(); }

    @Override
    public int position() { return mLeg.position(); }

    @Override
    public Integer startOffset() { return mLeg.startOffset(); }

    @Override
    public Integer endOffset() { return mLeg.endOffset(); }

    @Override
    public byte orientation() { return mLeg.orientation(); }

    @Override
    public String homology() { return mLeg.homology(); }

    @Override
    public Double alleleFrequency() { return mLeg.alleleFrequency(); }

    @Override
    public Integer inexactHomologyOffsetStart() { return mLeg.inexactHomologyOffsetStart(); }

    @Override
    public Integer inexactHomologyOffsetEnd() { return mLeg.inexactHomologyOffsetEnd(); }

    @Override
    public Integer tumorVariantFragmentCount() { return mLeg.tumorVariantFragmentCount(); }

    @Override
    public Integer tumorReferenceFragmentCount() { return mLeg.tumorReferenceFragmentCount(); }

    @Override
    public Integer normalVariantFragmentCount() { return mLeg.normalVariantFragmentCount(); }

    @Override
    public Integer normalReferenceFragmentCount() { return mLeg.normalReferenceFragmentCount(); }

    @Override
    public int anchoringSupportDistance() { return mLeg.anchoringSupportDistance(); }

    public double adjustedCopyNumber()
    {
        if(orientation() == 1)
        {
            return mLeftCopyNumber != null ? mLeftCopyNumber : 0;
        }
        else
        {
            return mRightCopyNumber != null ? mRightCopyNumber : 0;
        }
    }
}
