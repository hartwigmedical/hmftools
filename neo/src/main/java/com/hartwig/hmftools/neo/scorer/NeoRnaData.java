package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

public class NeoRnaData
{
    public final double TpmCancer;
    public final double TpmCohort;

    private double[] mTransExpression;
    private boolean mHasExpression;

    private int mFragmentSupport;
    private final int[] mBaseDepth;
    private boolean mHasCoverage;

    public NeoRnaData(double tpmCancer, double tpmCohort)
    {
        TpmCancer = tpmCancer;
        TpmCohort = tpmCohort;

        mFragmentSupport = 0;
        mBaseDepth = new int[] {0, 0};
        mTransExpression = new double[] {-1, -1}; // indicating not set
    }

    public boolean hasCoverage() { return mHasCoverage; }
    public boolean hasExpression() { return mHasExpression; }

    public double[] transExpression() { return mTransExpression; }
    public int fragmentSupport() { return mFragmentSupport; }
    public int[] baseDepth() { return mBaseDepth; }

    public void setExpression(final double[] expression)
    {
        mHasExpression = true;
        mTransExpression[FS_UP] = expression[FS_UP];
        mTransExpression[FS_DOWN] = expression[FS_DOWN];
    }

    public void setCoverage(int fragments, int depthUp, int depthDown)
    {
        mHasCoverage = true;
        mBaseDepth[FS_UP] = depthUp;
        mBaseDepth[FS_DOWN] = depthDown;
        mFragmentSupport = fragments;
    }

    public double getTPM()
    {
        if(mTransExpression[FS_UP] >= 0)
            return mTransExpression[FS_UP];

        return TpmCancer;
    }

    public double averageBaseDepth()
    {
        return (mBaseDepth[SE_START] + mBaseDepth[SE_END]) * 0.5;
    }
}
