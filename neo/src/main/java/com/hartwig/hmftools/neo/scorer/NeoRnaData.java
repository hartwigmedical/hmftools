package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

public class NeoRnaData
{
    public final double TpmCancer;
    public final double TpmCohort;
    public double[] TransExpression;
    public int FragmentSupport;
    public final int[] BaseDepth;

    public NeoRnaData(double tpmCancer, double tpmCohort)
    {
        TpmCancer = tpmCancer;
        TpmCohort = tpmCohort;
        FragmentSupport = 0;
        BaseDepth = new int[] {0, 0};
        TransExpression = new double[] {-1, -1}; // indicating not set
    }

    public double getTPM()
    {
        if(TransExpression[FS_UP] >= 0)
            return TransExpression[FS_UP];

        return TpmCancer;
    }

    public double baseDepth()
    {
        return (BaseDepth[SE_START] + BaseDepth[SE_END]) * 0.5;
    }
}
