package com.hartwig.hmftools.sage.sam;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface SamSlicerFactory {

    @NotNull
    SamSlicer create(@NotNull final GenomeRegion slice);

}
