package com.hartwig.hmftools.peach.panel;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;

import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class HaplotypePanel
{
    @NotNull
    private final Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel;

    public HaplotypePanel(@NotNull Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel)
    {
        this.geneToGeneHaplotypePanel = geneToGeneHaplotypePanel;
    }

    @NotNull
    public Map<Chromosome, Set<Integer>> getRelevantVariantPositions()
    {
        Map<Chromosome, Set<Integer>> chromosomeToRelevantPositions = new HashMap<>();
        for(GeneHaplotypePanel panel : geneToGeneHaplotypePanel.values())
        {
            for(Map.Entry<Chromosome, Set<Integer>> entry : panel.getRelevantVariantPositions().entrySet())
            {
                Chromosome chromosome = entry.getKey();
                Set<Integer> relevantPositions = entry.getValue();
                if(chromosomeToRelevantPositions.containsKey(chromosome))
                {
                    chromosomeToRelevantPositions.get(chromosome).addAll(relevantPositions);
                }
                else
                {
                    chromosomeToRelevantPositions.put(chromosome, new HashSet<>(relevantPositions));
                }
            }
        }
        return chromosomeToRelevantPositions;
    }

    public boolean isRelevantFor(@NotNull String eventId, String gene)
    {
        return isRelevantFor(HaplotypeEventFactory.fromId(eventId), gene);
    }

    public boolean isRelevantFor(@NotNull HaplotypeEvent event, String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).isRelevantFor(event);
    }

    @NotNull
    public List<NonDefaultHaplotype> getNonDefaultHaplotypes(@NotNull String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).nonDefaultHaplotypes;
    }

    @NotNull
    public DefaultHaplotype getDefaultHaplotype(@NotNull String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).defaultHaplotype;
    }

    @NotNull
    public String getWildTypeHaplotypeName(@NotNull String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).wildTypeHaplotypeName;
    }

    @NotNull
    public Set<String> getGenes()
    {
        return Sets.newHashSet(geneToGeneHaplotypePanel.keySet());
    }

    public int getHaplotypeCount()
    {
        return geneToGeneHaplotypePanel.values().stream().mapToInt(GeneHaplotypePanel::getHaplotypeCount).sum();
    }
}
