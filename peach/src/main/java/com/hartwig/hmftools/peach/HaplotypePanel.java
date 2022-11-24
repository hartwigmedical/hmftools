package com.hartwig.hmftools.peach;

import com.google.common.collect.Sets;
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

    public Set<String> getRelevantGenes(VariantHaplotypeEvent event)
    {
        return haplotypes.stream()
                .filter(h -> isRelevantHaplotype(h, event))
                .map(h -> h.gene)
                .collect(Collectors.toSet());
    }

    public Set<String> getGenes()
    {
        return haplotypes.stream().map(h -> h.gene).collect(Collectors.toSet());
    }

    private static boolean isRelevantHaplotype(Haplotype haplotype, VariantHaplotypeEvent event)
    {
        return haplotype.events.stream()
                .filter(e -> e instanceof VariantHaplotypeEvent)
                .map(e -> (VariantHaplotypeEvent) e)
                .map(VariantHaplotypeEvent::getCoveredPositions)
                .flatMap(Set::stream)
                .anyMatch(event.getCoveredPositions()::contains);
    }
}
