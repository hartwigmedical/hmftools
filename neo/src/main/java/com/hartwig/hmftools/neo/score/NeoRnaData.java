package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.INVALID_EXP;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

public class NeoRnaData
{
    private final double[] mTpmCancer;
    private final double[] mTpmPanCancer;

    private final double[] mTransExpression;
    private boolean mHasExpression;

    private int mFragmentSupport;
    private final int[] mBaseDepth;
    private boolean mHasCoverage;

    public static final double NO_TPM_VALUE = INVALID_EXP;

    public NeoRnaData()
    {
        mFragmentSupport = 0;
        mBaseDepth = new int[] {0, 0};
        mTransExpression = new double[] {NO_TPM_VALUE, NO_TPM_VALUE};
        mTpmCancer = new double[] {NO_TPM_VALUE, NO_TPM_VALUE};
        mTpmPanCancer = new double[] {NO_TPM_VALUE, NO_TPM_VALUE};
    }

    public boolean hasCoverage() { return mHasCoverage; }
    public boolean hasExpression() { return mHasExpression; }

    public double[] transExpression() { return mTransExpression; }
    public double[] tpmCancer() { return mTpmCancer; }
    public double[] tpmPanCancer() { return mTpmPanCancer; }

    public int fragmentSupport() { return mFragmentSupport; }
    public int[] baseDepth() { return mBaseDepth; }

    public void setCohortValues(final double[] tpmCancer, final double[] tpmPanCancer)
    {
        mTpmCancer[FS_UP] = tpmCancer[FS_UP];
        mTpmCancer[FS_DOWN] = tpmCancer[FS_DOWN];
        mTpmPanCancer[FS_UP] = tpmPanCancer[FS_UP];
        mTpmPanCancer[FS_DOWN] = tpmPanCancer[FS_DOWN];
    }

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

    public double getTPM(int stream)
    {
        if(mTransExpression[stream] != NO_TPM_VALUE)
            return mTransExpression[stream];

        return mTpmCancer[stream] != NO_TPM_VALUE ? mTpmCancer[stream] : mTpmPanCancer[stream];
    }

    public double averageBaseDepth()
    {
        return (mBaseDepth[SE_START] + mBaseDepth[SE_END]) * 0.5;
    }
}
