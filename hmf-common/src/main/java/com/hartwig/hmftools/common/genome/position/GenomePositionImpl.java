package com.hartwig.hmftools.common.genome.position;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public record GenomePositionImpl(String chromosome, int position) implements GenomePosition
{
    public GenomePositionImpl(HumanChromosome chromosome, RefGenomeVersion genomeVersion, int position)
    {
        this(genomeVersion.versionedChromosome(chromosome), position);
    }

    public GenomePositionImpl(GenomePosition position)
    {
        this(position.chromosome(), position.position());
    }
}
