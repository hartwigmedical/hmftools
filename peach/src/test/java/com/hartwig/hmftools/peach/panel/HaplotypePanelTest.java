package com.hartwig.hmftools.peach.panel;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;

import org.junit.Test;

public class HaplotypePanelTest
{
    @Test
    public void testGetRelevantVariantPositionsSameChromosome()
    {
        ImmutableList<HaplotypeEvent> eventsToIgnore = ImmutableList.copyOf(Collections.emptyList());

        VariantHaplotypeEvent fake1HaplotypeEvent = new VariantHaplotypeEvent(HumanChromosome._1, 1000, "A", "T");
        NonDefaultHaplotype fake1NonDefaultHaplotype =
                new NonDefaultHaplotype("*2", false, ImmutableList.copyOf(List.of(fake1HaplotypeEvent)));
        DefaultHaplotype fake1DefaultHaplotype = new DefaultHaplotype("*1", true, eventsToIgnore);
        GeneHaplotypePanel fake1Panel =
                new GeneHaplotypePanel(fake1DefaultHaplotype, ImmutableList.copyOf(List.of(fake1NonDefaultHaplotype)), "*1");

        VariantHaplotypeEvent fake2HaplotypeEvent = new VariantHaplotypeEvent(HumanChromosome._1, 3000, "AG", "A");
        NonDefaultHaplotype fake2NonDefaultHaplotype =
                new NonDefaultHaplotype("*2", false, ImmutableList.copyOf(List.of(fake2HaplotypeEvent)));
        DefaultHaplotype fake2DefaultHaplotype = new DefaultHaplotype("*1", true, eventsToIgnore);
        GeneHaplotypePanel fake2Panel =
                new GeneHaplotypePanel(fake2DefaultHaplotype, ImmutableList.copyOf(List.of(fake2NonDefaultHaplotype)), "*1");
        HaplotypePanel panel = new HaplotypePanel(Map.of("FAKE1", fake1Panel, "FAKE2", fake2Panel));

        Map<Chromosome, Set<Integer>> relevantVariantPositionMap = panel.getRelevantVariantPositions();

        assertEquals(1, relevantVariantPositionMap.keySet().size());
        assertTrue(relevantVariantPositionMap.containsKey(HumanChromosome._1));

        Set<Integer> relevantPositions = relevantVariantPositionMap.get(HumanChromosome._1);
        assertEquals(3, relevantPositions.size());
        assertTrue(relevantPositions.contains(1000));
        assertTrue(relevantPositions.contains(3000));
        assertTrue(relevantPositions.contains(3001));
    }
}
