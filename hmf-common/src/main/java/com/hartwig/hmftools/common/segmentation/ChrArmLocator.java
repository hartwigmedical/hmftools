package com.hartwig.hmftools.common.segmentation;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public interface ChrArmLocator
{
    static ChrArmLocator defaultLocator(RefGenomeVersion refGenomeVersion)
    {
        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenomeVersion);
        return cobaltRatio ->
        {
            int centre = refGenomeCoordinates.centromere(cobaltRatio.chromosome());
            HumanChromosome chromosome = HumanChromosome.fromString(cobaltRatio.chromosome());
            Arm arm = cobaltRatio.position() < centre ? Arm.P : Arm.Q;
            return new ChrArm(chromosome, arm);
        };
    }

    ChrArm map(GenomePosition cobaltRatio);
}
