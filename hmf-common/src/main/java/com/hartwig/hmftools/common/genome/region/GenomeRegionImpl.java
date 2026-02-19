package com.hartwig.hmftools.common.genome.region;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public record GenomeRegionImpl(String chromosome, int start, int end) implements GenomeRegion
{
    public GenomeRegionImpl(HumanChromosome humanChromosome, RefGenomeVersion genomeVersion, int start, int end)
    {
        this(genomeVersion.versionedChromosome(humanChromosome), start, end);
    }
}