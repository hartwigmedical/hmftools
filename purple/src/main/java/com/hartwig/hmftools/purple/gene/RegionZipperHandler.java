package com.hartwig.hmftools.purple.gene;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface RegionZipperHandler<S extends GenomeRegion, T extends GenomeRegion>
{
    void enterChromosome(@NotNull String chromosome);

    void primary(@NotNull S region);

    void secondary(@NotNull T region);
}