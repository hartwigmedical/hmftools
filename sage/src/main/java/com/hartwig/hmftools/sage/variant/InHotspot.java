package com.hartwig.hmftools.sage.variant;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

class InHotspot {

    private static final int DISTANCE = 5;
    private final Multimap<Chromosome, VariantHotspot> hotspots;

    InHotspot(@NotNull final Chromosome chromosome,@NotNull final List<VariantHotspot> hotspots) {
        this.hotspots = ArrayListMultimap.create();
        this.hotspots.putAll(chromosome, hotspots);
    }

    InHotspot(@NotNull final Multimap<Chromosome, VariantHotspot> hotspots) {
        this.hotspots = hotspots;
    }

    public boolean isOnHotspot(@NotNull final VariantHotspot context) {
        if (HumanChromosome.contains(context.chromosome())) {
            final Chromosome chromosome = HumanChromosome.fromString(context.chromosome());
            if (hotspots.containsKey(chromosome)) {
                return hotspots.get(chromosome).stream().anyMatch(x -> exactMatch(x, context));
            }
        }

        return false;
    }

    public boolean isNearHotspot(@NotNull final VariantHotspot context) {
        if (HumanChromosome.contains(context.chromosome())) {
            final Chromosome chromosome = HumanChromosome.fromString(context.chromosome());
            if (hotspots.containsKey(chromosome)) {
                return hotspots.get(chromosome).stream().anyMatch(x -> overlaps(x, context));
            }
        }

        return false;
    }

    private static boolean overlaps(@NotNull final VariantHotspot hotspot, @NotNull final VariantHotspot target) {
        long variantStart = target.position();
        long variantEnd = target.position() + target.ref().length() - 1 + DISTANCE;

        long ponStart = hotspot.position();
        long ponEnd = hotspot.position() + hotspot.ref().length() - 1 + DISTANCE;

        return variantStart <= ponEnd && variantEnd >= ponStart;
    }

    private static boolean exactMatch(@NotNull final VariantHotspot hotspot, @NotNull final VariantHotspot target) {
        return hotspot.position() == target.position() && hotspot.ref().equals(target.ref()) && hotspot.alt().equals(target.alt());
    }
}

