package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class HaplotypePanel
{
    @NotNull
    private final List<Haplotype> haplotypes;

    public HaplotypePanel(@NotNull List<Haplotype> haplotypes)
    {
        this.haplotypes = List.copyOf(haplotypes);
    }

    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        Stream<VariantHaplotypeEvent> variantHaplotypeEvents = haplotypes.stream()
                .filter(h -> ! h.wildType) // Ignore variants configured for wild type haplotypes, since they should be ignored anyway
                .map(h -> h.events)
                .flatMap(Collection::stream)
                .filter(e -> e instanceof VariantHaplotypeEvent)
                .map(VariantHaplotypeEvent.class::cast);
        return variantHaplotypeEvents.collect(
                Collectors.groupingBy(
                        e -> e.chromosome,
                        Collectors.mapping(VariantHaplotypeEvent::getCoveredPositions, Collectors.flatMapping(Collection::stream, Collectors.toSet()))
                )
        );
    }
}
