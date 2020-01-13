package com.hartwig.hmftools.sage.sam;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

public class SamSlicerFactoryMultiImpl implements SamSlicerFactory {

    private final SageConfig config;
    private final ListMultimap<Chromosome, GenomeRegion> panel;

    public SamSlicerFactoryMultiImpl(@NotNull final SageConfig config, @NotNull final ListMultimap<Chromosome, GenomeRegion> panel) {
        this.config = config;
        this.panel = panel;
    }

    @NotNull
    public SamSlicer create(@NotNull final GenomeRegion slice) {
        return config.panelOnly() ? panelOnly(slice) : fullSlice(slice);
    }

    @NotNull
    private SamSlicer fullSlice(@NotNull final GenomeRegion slice) {
        return new SamSlicer(config.minMapQuality(), slice);
    }

    @NotNull
    private SamSlicer panelOnly(@NotNull final GenomeRegion slice) {
        Chromosome chromosome = HumanChromosome.fromString(slice.chromosome());
        return new SamSlicer(config.minMapQuality(), slice, panel.get(chromosome));
    }

}
