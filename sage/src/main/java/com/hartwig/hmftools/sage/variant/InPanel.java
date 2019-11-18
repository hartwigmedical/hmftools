package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.select.RegionSelector;

import org.jetbrains.annotations.NotNull;

class InPanel {

    private final Map<Chromosome, RegionSelector<GenomeRegion>> panelSelectorMap;

    InPanel(@NotNull final Chromosome chromosome, @NotNull final List<GenomeRegion> regions) {
        this.panelSelectorMap = Maps.newHashMap();
        panelSelectorMap.put(chromosome, new RegionSelector<>(regions));
    }

    InPanel(@NotNull final ListMultimap<Chromosome, GenomeRegion> panelRegions) {
        this.panelSelectorMap = Maps.newHashMap();
        for (Chromosome chromosome : panelRegions.keySet()) {
            List<GenomeRegion> regions = panelRegions.get(chromosome);
            panelSelectorMap.put(chromosome, new RegionSelector<>(regions));
        }
    }

    public boolean inPanel(@NotNull final GenomePosition position) {
        return HumanChromosome.contains(position.chromosome()) && Optional.ofNullable(panelSelectorMap.get(HumanChromosome.fromString(
                position.chromosome()))).flatMap(x -> x.select(position.position())).isPresent();

    }

}
