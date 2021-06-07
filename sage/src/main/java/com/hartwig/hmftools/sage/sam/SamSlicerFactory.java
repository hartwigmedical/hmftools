package com.hartwig.hmftools.sage.sam;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

public class SamSlicerFactory
{

    private final SageConfig config;
    private final List<GenomeRegion> panel;

    public SamSlicerFactory(@NotNull final SageConfig config, @NotNull final List<GenomeRegion> panel)
    {
        this.config = config;
        this.panel = panel;
    }

    @NotNull
    public SamSlicer create(@NotNull final GenomeRegion slice)
    {
        return config.panelOnly() ? panelOnly(slice) : fullSlice(slice);
    }

    @NotNull
    private SamSlicer fullSlice(@NotNull final GenomeRegion slice)
    {
        return new SamSlicer(0, slice);
    }

    @NotNull
    private SamSlicer panelOnly(@NotNull final GenomeRegion slice)
    {
        return new SamSlicer(0, slice, panel);
    }
}
