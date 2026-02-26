package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;

public record ChrArmCopyNumber(
        HumanChromosome chromosome, Arm arm, double meanCopyNumber, double medianCopyNumber, double minCopyNumber, double maxCopyNumber)
{
    public boolean includeInReport()
    {
        if(chromosome.hasShortArm())
            return arm == Arm.Q;

        return true;
    }
}
