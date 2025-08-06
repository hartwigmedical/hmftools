package com.hartwig.hmftools.cobalt.e2e;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public record GenomeBackedGcFileSection(HumanChromosome chromosome, int start, int stop, RefGenomeInterface genome) implements GcFileSection
{
    @Override
    public double gcRatio(final int position)
    {
        final String chrName = RefGenomeVersion.V38.versionedChromosome(chromosome.toString());
        String bases = genome.getBaseString(chrName, position, position + 1000).toUpperCase(); // TODO not always 1000
        int gcCount = 0;
        for (char chr : bases.toCharArray()) {
            if (chr == 'G') gcCount++;
            if (chr == 'C') gcCount++;
        }
        return gcCount / (double) bases.length();
    }
}
