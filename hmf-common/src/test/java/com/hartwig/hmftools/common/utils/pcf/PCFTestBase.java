package com.hartwig.hmftools.common.utils.pcf;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PCFTestBase
{
    ChrBaseRegion region(HumanChromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(V38.versionedChromosome(chromosome), start, end);
    }
}
