package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GeneHaplotypePanel
{
    @NotNull
    private final String gene;
    @NotNull
    private final Haplotype wildTypeHaplotype;
    @NotNull
    private final List<Haplotype> nonWildTypeHaplotypes;

    public GeneHaplotypePanel(
            @NotNull String gene,
            @NotNull Haplotype wildTypeHaplotype,
            @NotNull List<Haplotype> nonWildTypeHaplotypes
    )
    {
        this.gene = gene;
        this.wildTypeHaplotype = wildTypeHaplotype;
        this.nonWildTypeHaplotypes = List.copyOf(nonWildTypeHaplotypes);
    }

    public static GeneHaplotypePanel fromHaplotypes(@NotNull String gene, @NotNull List<Haplotype> haplotypes)
    {
        List<Haplotype> wildTypeHaplotypes = haplotypes.stream().filter(h -> h.wildType).collect(Collectors.toList());
        if (wildTypeHaplotypes.size() != 1)
        {
            throw new RuntimeException(String.format("Gene '%s' does not have precisely one wild type haplotype", gene));
        }
        Haplotype wildTypeHaplotype = wildTypeHaplotypes.get(0);
        List<Haplotype> nonWildTypeHaplotypes = haplotypes.stream().filter(h -> !h.wildType).collect(Collectors.toList());
        return new GeneHaplotypePanel(gene, wildTypeHaplotype, nonWildTypeHaplotypes);
    }

    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        Stream<VariantHaplotypeEvent> variantHaplotypeEvents = nonWildTypeHaplotypes.stream()
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
        return nonWildTypeHaplotypes.stream().anyMatch(h -> h.isRelevantFor(event));
    }
}
