package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.FittingRegion;

record FittedRegionImpl(HumanChromosome chr, int start, int end,
                        GermlineStatus germlineStatus, int bafCount, double observedBAF,
                        double observedNormalRatio, double observedTumorRatio) implements FittingRegion
{
    @Override
    public String chromosome()
    {
        return V38.versionedChromosome(chr);
    }
}
