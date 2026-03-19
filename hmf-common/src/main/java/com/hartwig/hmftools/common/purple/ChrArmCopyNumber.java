package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.driver.DriverInterpretation;
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

    public static final double HIGH_DRIVER_THRESHOLD = 0.4;
    public static final double LOW_DRIVER_THRESHOLD = 0.25;

    public static final double HIGH_DRIVER_ARM_GAIN_THRESHOLD = 1 + HIGH_DRIVER_THRESHOLD;
    public static final double LOW_DRIVER_ARM_CN_GAIN_THRESHOLD = 1 + LOW_DRIVER_THRESHOLD;
    public static final double HIGH_DRIVER_ARM_LOSS_THRESHOLD = 1 - HIGH_DRIVER_THRESHOLD;
    public static final double LOW_DRIVER_ARM_CN_LOSS_THRESHOLD = 1 - LOW_DRIVER_THRESHOLD;

    public static DriverInterpretation driverInterpretation(final ChrArmCopyNumber chrArmCopyNumber, final double ploidy)
    {
        if((chrArmCopyNumber.meanCopyNumber() > HIGH_DRIVER_ARM_GAIN_THRESHOLD * ploidy)
        || (chrArmCopyNumber.meanCopyNumber() < HIGH_DRIVER_ARM_LOSS_THRESHOLD * ploidy))
        {
            return DriverInterpretation.HIGH;
        }
        if((chrArmCopyNumber.meanCopyNumber() > LOW_DRIVER_ARM_CN_GAIN_THRESHOLD * ploidy)
        || (chrArmCopyNumber.meanCopyNumber() < LOW_DRIVER_ARM_CN_LOSS_THRESHOLD * ploidy))
        {
            return DriverInterpretation.LOW;
        }

        return null;
    }
}
