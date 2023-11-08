package com.hartwig.hmftools.peach.panel;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GeneHaplotypePanel
{
    @NotNull
    public final DefaultHaplotype defaultHaplotype;
    @NotNull
    public final ImmutableList<NonDefaultHaplotype> nonDefaultHaplotypes;
    @NotNull
    public final String wildTypeHaplotypeName;

    public GeneHaplotypePanel(
            @NotNull DefaultHaplotype defaultHaplotype,
            @NotNull ImmutableList<NonDefaultHaplotype> nonDefaultHaplotypes,
            @NotNull String wildTypeHaplotypeName
    )
    {
        this.defaultHaplotype = defaultHaplotype;
        this.nonDefaultHaplotypes = nonDefaultHaplotypes;
        this.wildTypeHaplotypeName = wildTypeHaplotypeName;
    }

    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        Stream<VariantHaplotypeEvent> variantHaplotypeEvents = nonDefaultHaplotypes.stream()
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

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        boolean isEventToIgnore = defaultHaplotype.eventsToIgnore.stream()
                .map(HaplotypeEvent::id)
                .anyMatch(e -> e.equals(event.id()));
        return !isEventToIgnore && nonDefaultHaplotypes.stream().anyMatch(h -> h.isRelevantFor(event));
    }

    public int getHaplotypeCount()
    {
        return nonDefaultHaplotypes.size() + 1;
    }
}
