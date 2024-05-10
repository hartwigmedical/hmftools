package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PanelLoaderTest
{
    @Test
    public void testLoad()
    {
        String filePath = getTestResourcePath("haplotypes.complicated.37.tsv");
        HaplotypePanel panel = PanelLoader.loadHaplotypePanel(filePath);

        assertEquals(2, panel.getGenes().size());

        assertEqualDefaultHaplotypes(new DefaultHaplotype("*9A", false, ImmutableList.of()), panel.getDefaultHaplotype("DPYD"));
        assertEqualDefaultHaplotypes(new DefaultHaplotype("*1", true, ImmutableList.of()), panel.getDefaultHaplotype("UGT1A1"));

        assertEquals(4, panel.getNonDefaultHaplotypes("UGT1A1").size());
        NonDefaultHaplotype expectedFirstUgt1a1NonDefaultHaplotype = new NonDefaultHaplotype("*28", false,
                ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._2, 234668879, "C", "CAT"))
        );
        assertEqualNonDefaultHaplotypes(expectedFirstUgt1a1NonDefaultHaplotype, panel.getNonDefaultHaplotypes("UGT1A1").get(0));

        assertEquals(6, panel.getNonDefaultHaplotypes("DPYD").size());
        NonDefaultHaplotype expectedSecondDpydNonDefaultHaplotype = new NonDefaultHaplotype("*2A", false,
            ImmutableList.of(
                    new VariantHaplotypeEvent(HumanChromosome._1, 97915614, "C", "T"),
                    new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
            )
        );
        assertEqualNonDefaultHaplotypes(expectedSecondDpydNonDefaultHaplotype, panel.getNonDefaultHaplotypes("DPYD").get(1));
    }

    private static void assertEqualDefaultHaplotypes(@NotNull DefaultHaplotype expected, @NotNull DefaultHaplotype actual)
    {
        assertEquals(expected.getName(), actual.getName());
        assertEquals(expected.isWildType(), actual.isWildType());
        assertEqualHaplotypeEventLists(expected.eventsToIgnore, actual.eventsToIgnore);
    }

    private static void assertEqualNonDefaultHaplotypes(@NotNull NonDefaultHaplotype expected, @NotNull NonDefaultHaplotype actual)
    {
        assertEquals(expected.getName(), actual.getName());
        assertEquals(expected.isWildType(), actual.isWildType());
        assertEqualHaplotypeEventLists(expected.events, actual.events);
    }

    private static void assertEqualHaplotypeEventLists(@NotNull ImmutableList<HaplotypeEvent> expected, @NotNull ImmutableList<HaplotypeEvent> actual)
    {
        assertEquals(expected.size(), actual.size());
        for(int i = 0; i < expected.size(); i++)
        {
            assertEqualHaplotypeEvents(expected.get(i), actual.get(i));
        }
    }

    private static void assertEqualHaplotypeEvents(@NotNull HaplotypeEvent expected, @NotNull HaplotypeEvent actual)
    {
        assertEquals(expected.id(), actual.id());
    }
}
