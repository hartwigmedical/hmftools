package com.hartwig.hmftools.peach.panel;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class HaplotypePanel
{
    private final Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel;

    public HaplotypePanel(final Map<String, GeneHaplotypePanel> geneToGeneHaplotypePanel)
    {
        this.geneToGeneHaplotypePanel = geneToGeneHaplotypePanel;
    }

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

    public boolean isRelevantFor(final String eventId, String gene)
    {
        return isRelevantFor(HaplotypeEventFactory.fromId(eventId), gene);
    }

    public boolean isRelevantFor(final HaplotypeEvent event, String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).isRelevantFor(event);
    }

    public List<NonDefaultHaplotype> getNonDefaultHaplotypes(final String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).nonDefaultHaplotypes;
    }

    public DefaultHaplotype getDefaultHaplotype(final String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).defaultHaplotype;
    }

    public String getWildTypeHaplotypeName(final String gene)
    {
        return geneToGeneHaplotypePanel.get(gene).wildTypeHaplotypeName;
    }

    public Set<String> getGenes()
    {
        return Sets.newHashSet(geneToGeneHaplotypePanel.keySet());
    }

    public int getHaplotypeCount()
    {
        return geneToGeneHaplotypePanel.values().stream().mapToInt(GeneHaplotypePanel::getHaplotypeCount).sum();
    }
}
