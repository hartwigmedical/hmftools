package com.hartwig.hmftools.linx.cn;

import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

public class CnArmStats
{
    public final String Chromosome;
    public final ChromosomeArm Arm;

    public int SegmentCount;
    public double AverageCopyNumber;
    public double MedianCopyNumber;
    public double MaxCopyNumber;
    public double MinCopyNumber;
    public boolean HasLOH;
    public double TelomereCopyNumber;
    public double CentromereCopyNumber;

    public CnArmStats(final String chromosome, final ChromosomeArm arm)
    {
        Chromosome = chromosome;
        Arm = arm;

        SegmentCount = 0;
        AverageCopyNumber = 0;
        MedianCopyNumber = 0;
        MaxCopyNumber = 0;
        MinCopyNumber = -1;
        HasLOH = false;
        TelomereCopyNumber = 0;
        CentromereCopyNumber = 0;
    }
}
