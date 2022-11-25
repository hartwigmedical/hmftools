package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
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

    public Set<String> getRelevantGenes(String eventId)
    {
        return geneToGeneHaplotypePanel.entrySet().stream()
                .filter(e -> e.getValue().isRelevantFor(eventId))
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());
    }

    public boolean isRelevantFor(String eventId, String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).isRelevantFor(eventId);
    }
}
