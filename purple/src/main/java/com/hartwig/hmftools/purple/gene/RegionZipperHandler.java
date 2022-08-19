package com.hartwig.hmftools.purple.gene;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface RegionZipperHandler<S extends GenomeRegion, T extends GenomeRegion>
{
    void primary(final S region);

    void secondary(final T region);
}