package com.hartwig.hmftools.cobalt.segmentation;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

interface ChrArmLocator
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

    ChrArm map(CobaltRatio cobaltRatio);
}
