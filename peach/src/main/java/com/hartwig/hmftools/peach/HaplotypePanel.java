package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Collectors;

public class HaplotypePanel
{
    @NotNull
    private final Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel;

    public HaplotypePanel(@NotNull Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel)
    {
        this.geneToGeneHaplotypePanel = geneToGeneHaplotypePanel;
    }

    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        return geneToGeneHaplotypePanel.values().stream()
                .map(GeneHaplotypePanel::getRelevantVariantPositions)
                .flatMap(m -> m.entrySet().stream())
                .collect(
                        Collectors.groupingBy(
                                Map.Entry::getKey,
                                Collectors.mapping(Map.Entry::getValue, Collectors.flatMapping(Collection::stream, Collectors.toSet()))
                        )
                );
    }

    public Set<String> getRelevantGenes(HaplotypeEvent event)
    {
        return geneToGeneHaplotypePanel.entrySet().stream()
                .filter(e -> e.getValue().isRelevantFor(event))
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());
    }

    public Set<String> getGenes()
    {
        return new HashSet<>(geneToGeneHaplotypePanel.keySet());
    }


}
