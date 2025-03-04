package com.hartwig.hmftools.peach.panel;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class GeneHaplotypePanel
{
    public final DefaultHaplotype defaultHaplotype;
    public final ImmutableList<NonDefaultHaplotype> nonDefaultHaplotypes;
    public final String wildTypeHaplotypeName;

    public GeneHaplotypePanel(
            final DefaultHaplotype defaultHaplotype, final ImmutableList<NonDefaultHaplotype> nonDefaultHaplotypes,
            final String wildTypeHaplotypeName)
    {
        this.defaultHaplotype = defaultHaplotype;
        this.nonDefaultHaplotypes = nonDefaultHaplotypes;
        this.wildTypeHaplotypeName = wildTypeHaplotypeName;
    }

    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        Map<Chromosome, Set<Integer>> chromosomeToRelevantVariantPositions = new HashMap<>();
        for(VariantHaplotypeEvent event : getNonDefaultVariantHaplotypeEvents())
        {
            Chromosome chromosome = event.chromosome;
            Set<Integer> relevantPositions = event.getCoveredPositions();
            if(chromosomeToRelevantVariantPositions.containsKey(chromosome))
            {
                chromosomeToRelevantVariantPositions.get(chromosome).addAll(relevantPositions);
            }
            else
            {
                chromosomeToRelevantVariantPositions.put(chromosome, new HashSet<>(relevantPositions));
            }
        }
        return chromosomeToRelevantVariantPositions;
    }

    public boolean isRelevantFor(final HaplotypeEvent event)
    {
        boolean isEventToIgnore = defaultHaplotype.eventsToIgnore.stream().map(HaplotypeEvent::id).anyMatch(e -> e.equals(event.id()));
        return !isEventToIgnore && nonDefaultHaplotypes.stream().anyMatch(h -> h.isRelevantFor(event));
    }

    public int getHaplotypeCount()
    {
        return nonDefaultHaplotypes.size() + 1;
    }

    private List<VariantHaplotypeEvent> getNonDefaultVariantHaplotypeEvents()
    {
        return nonDefaultHaplotypes.stream()
                .map(h -> h.events)
                .flatMap(Collection::stream)
                .filter(e1 -> e1 instanceof VariantHaplotypeEvent)
                .map(VariantHaplotypeEvent.class::cast)
                .collect(Collectors.toList());
    }
}
