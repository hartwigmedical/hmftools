package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;
import static com.hartwig.hmftools.peach.panel.TestHaplotypePanelFactory.createDefaultTestHaplotypePanel;
import static com.hartwig.hmftools.peach.panel.TestHaplotypePanelFactory.createDpydTestHaplotypePanel;
import static com.hartwig.hmftools.peach.panel.TestHaplotypePanelFactory.createUgt1a1TestHaplotypePanel;

import static junit.framework.TestCase.assertEquals;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.panel.GeneHaplotypePanel;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PeachAlgoTest
{
    @Test
    public void testNoPanelOrEvents()
    {
        HaplotypePanel haplotypePanel = new HaplotypePanel(new HashMap<>());
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, new HashMap<>());
        assertEquals(0, geneToHaplotypeAnalysis.size());
    }

    @Test
    public void testNoEvents()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, new HashMap<>());

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*9A", 2))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*9A", 2))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testHomozygousRefEvents()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr1_98348885_G_A", 0,  // relevant event
                "VAR_chr1_98348885_G_T", 0 // overlapping irrelevant event
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                eventIdToCount,
                List.of(new HaplotypeCombination(Map.of("*9A", 2))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*9A", 2))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testNoCallRelevantEvent()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = new HashMap<>();
        eventIdToCount.put("VAR_chr1_98348885_G_T", null);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                eventIdToCount,
                new ArrayList<>(),
                "*9A",
                "*1",
                PeachQCStatus.FAIL_EVENT_WITH_UNKNOWN_COUNT,
                null
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testNoCallOverlappingEvent()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = new HashMap<>();
        eventIdToCount.put("VAR_chr1_98348885_G_T", null);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                eventIdToCount,
                new ArrayList<>(),
                "*9A",
                "*1",
                PeachQCStatus.FAIL_EVENT_WITH_UNKNOWN_COUNT,
                null
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testWildTypes()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of("VAR_chr1_98348885_G_A", 2);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("VAR_chr1_98348885_G_A", 2),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testHeterozygousNonWildType()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr1_98348885_G_A", 2,
                "VAR_chr1_98205966_GATGA_G", 1,
                "VAR_chr1_98039419_C_T", 1,
                "VAR_chr1_98045449_G_C", 1,
                "VAR_chr2_234668879_C_CAT", 1,
                "VAR_chr2_234668879_CAT_C", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr1_98348885_G_A", 2,
                        "VAR_chr1_98205966_GATGA_G", 1,
                        "VAR_chr1_98039419_C_T", 1,
                        "VAR_chr1_98045449_G_C", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*7", 1, "*B3", 1))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*7", 1, "*B3", 1))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CAT", 1,
                        "VAR_chr2_234668879_CAT_C", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*28", 1, "*36", 1))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*28", 1, "*36", 1))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testHeterozygousWildType()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr1_98348885_G_A", 1,
                "VAR_chr2_234669144_G_A", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("VAR_chr1_98348885_G_A", 1),
                List.of(new HaplotypeCombination(Map.of("*9A", 1, "*1", 1))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*9A", 1, "*1", 1))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("VAR_chr2_234669144_G_A", 1),
                List.of(new HaplotypeCombination(Map.of("*1", 1, "*6", 1))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 1, "*6", 1))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testHeterozygousDefaultHaplotype()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr1_98348885_G_A", 1,
                "VAR_chr1_97915614_C_T", 1,
                "VAR_chr2_234668879_C_CATAT", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr1_98348885_G_A", 1,
                        "VAR_chr1_97915614_C_T", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*2A", 1, "*9A", 1))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*2A", 1, "*9A", 1))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CATAT", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*1", 1, "*37", 1))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 1, "*37", 1))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testHomozygousNonWildType()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr1_97981343_A_C", 2,
                "VAR_chr1_98348885_G_A", 2,
                "VAR_chr2_234668879_C_CAT", 2
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr1_98348885_G_A", 2,
                        "VAR_chr1_97981343_A_C", 2
                ),
                List.of(new HaplotypeCombination(Map.of("*13", 2))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*13", 2))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CAT", 2
                ),
                List.of(new HaplotypeCombination(Map.of("*28", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*28", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testTooManyAllelesUgt1a1()
    {
        HaplotypePanel haplotypePanel = createDefaultTestHaplotypePanel();
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr2_234668879_C_CAT", 2,
                "VAR_chr2_234669144_G_A", 2
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(),
                List.of(new HaplotypeCombination(Map.of("*9A", 2))),
                "*9A",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*9A", 2))
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CAT", 2,
                        "VAR_chr2_234669144_G_A", 2
                ),
                List.of(new HaplotypeCombination(Map.of("*28", 2, "*6", 2))),
                "*1",
                "*1",
                PeachQCStatus.WARN_TOO_MANY_ALLELES_FOUND,
                new HaplotypeCombination(Map.of("*28", 2, "*6", 2))
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));

        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }

    @Test
    public void testOverlappingGenes()
    {
        HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(getTestResourcePath("haplotypes.overlap.tsv"));
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr3_1234567_CAT_C", 1,
                "VAR_chr3_1234567_C_CAT", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedFake1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_CAT_C", 1,
                        "VAR_chr3_1234567_C_CAT", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*2", 1, "*4", 1))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*2", 1, "*4", 1))
        );
        HaplotypeAnalysis expectedFake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_C_CAT", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*2", 1, "*1", 1))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*2", 1, "*1", 1))
        );
        HaplotypeAnalysis expectedFake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_CAT_C", 1,
                        "VAR_chr3_1234567_C_CAT", 1
                ),
                new ArrayList<>(),
                "*1",
                "*1",
                PeachQCStatus.FAIL_NO_COMBINATION_FOUND,
                null
        );
        assertEquals(3, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE2"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE3"));

        assertEqualHaplotypeAnalysis(expectedFake1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
        assertEqualHaplotypeAnalysis(expectedFake2HaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE2"));
        assertEqualHaplotypeAnalysis(expectedFake3HaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE3"));
    }

    @Test
    public void testAmbiguousCallingPreferFewerNonWildTypeHaplotypes()
    {
        HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(getTestResourcePath("haplotypes.ambiguous.tsv"));
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr3_1234567_C_CAT", 2,
                "VAR_chr3_1234765_G_C", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_C_CAT", 2,
                        "VAR_chr3_1234765_G_C", 1
                ),
                List.of(
                        new HaplotypeCombination(Map.of("*2", 1, "*3", 1)),
                        new HaplotypeCombination(Map.of("*3", 2, "*4", 1))
                ),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*2", 1, "*3", 1))
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
    }

    @Test
    public void testAmbiguousCallingPreferFewerHaplotypes()
    {
        // strictly prefer two non-wild type haplotypes over two wild type and two non-wild type haplotypes.
        HaplotypeEvent event1 = new VariantHaplotypeEvent(HumanChromosome._1, 100, "A", "C");
        HaplotypeEvent event2 = new VariantHaplotypeEvent(HumanChromosome._1, 200, "A", "ATT");
        DefaultHaplotype defaultHaplotype = new DefaultHaplotype("*2", false, ImmutableList.of());
        ImmutableList<NonDefaultHaplotype> nonDefaultHaplotypes = ImmutableList.of(
                new NonDefaultHaplotype("*1", true, ImmutableList.of(event1)),
                new NonDefaultHaplotype("*3", true, ImmutableList.of(event2)),
                new NonDefaultHaplotype("*4", true, ImmutableList.of(event1, event2))
        );
        GeneHaplotypePanel geneHaplotypePanel = new GeneHaplotypePanel(defaultHaplotype, nonDefaultHaplotypes, "*1");
        HaplotypePanel haplotypePanel = new HaplotypePanel(Map.of("FAKE1", geneHaplotypePanel));

        Map<String, Integer> eventIdToCount = Map.of(
                event1.id(), 2,
                event2.id(), 2
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr1_100_A_C", 2,
                        "VAR_chr1_200_A_ATT", 2
                ),
                List.of(
                        new HaplotypeCombination(Map.of("*1", 2, "*3", 2)),
                        new HaplotypeCombination(Map.of("*1", 1, "*3", 1, "*4", 1)),
                        new HaplotypeCombination(Map.of("*4", 2))
                ),
                "*2",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*4", 2))
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
    }

    @Test
    public void testAmbiguousCallingNoBestHaplotype()
    {
        HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(getTestResourcePath("haplotypes.ambiguous.tsv"));
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr3_1234567_C_CAT", 1,
                "VAR_chr3_1234765_G_C", 1,
                "VAR_chr3_1234777_AA_T", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_C_CAT", 1,
                        "VAR_chr3_1234765_G_C", 1,
                        "VAR_chr3_1234777_AA_T", 1
                ),
                List.of(
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1)),
                        new HaplotypeCombination(Map.of("*3", 1, "*6", 1)),
                        new HaplotypeCombination(Map.of("*3", 1, "*4", 1, "*5", 1))
                ),
                "*1",
                "*1",
                PeachQCStatus.FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND,
                null
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
    }

    @Test
    public void testAmbiguousCallingPreferWildType()
    {
        HaplotypePanel haplotypePanel = PanelLoader.loadHaplotypePanel(getTestResourcePath("haplotypes.ambiguous.tsv"));
        Map<String, Integer> eventIdToCount = Map.of(
                "VAR_chr3_1234567_C_CAT", 1,
                "VAR_chr3_1234765_G_C", 1
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(haplotypePanel, eventIdToCount);

        HaplotypeAnalysis expectedHaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_C_CAT", 1,
                        "VAR_chr3_1234765_G_C", 1
                ),
                List.of(
                        new HaplotypeCombination(Map.of("*2", 1, "*1", 1)),
                        new HaplotypeCombination(Map.of("*3", 1, "*4", 1))
                ),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*2", 1, "*1", 1))
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
    }

    @Test
    public void testSimpleHomozygousDpydHaplotypeCombinations()
    {
        testAllSimpleHomozygousHaplotypeCombinations("DPYD", createDpydTestHaplotypePanel());
    }

    @Test
    public void testSimpleHeterozygousDpydHaplotypeCombinations()
    {
        testAllSimpleHeterozygousHaplotypeCombinations("DPYD", createDpydTestHaplotypePanel());
    }

    @Test
    public void testSimpleHomozygousUgt1a1HaplotypeCombinations()
    {
        testAllSimpleHomozygousHaplotypeCombinations("UGT1A1", createUgt1a1TestHaplotypePanel());
    }

    @Test
    public void testSimpleHeterozygousUgt1a1HaplotypeCombinations()
    {
        testAllSimpleHeterozygousHaplotypeCombinations("UGT1A1", createUgt1a1TestHaplotypePanel());
    }

    private static void testAllSimpleHomozygousHaplotypeCombinations(String gene, HaplotypePanel panel)
    {
        Map<String, List<String>> haplotypeNameToEventIds = panel.getNonDefaultHaplotypes(gene)
                .stream()
                .collect(Collectors.toMap(NonDefaultHaplotype::getName, PeachAlgoTest::getEventIds));
        String defaultHaplotypeName = panel.getDefaultHaplotype(gene).getName();
        haplotypeNameToEventIds.put(defaultHaplotypeName, Collections.emptyList());

        List<String> haplotypeNames = haplotypeNameToEventIds.keySet().stream().sorted().collect(Collectors.toList());
        for(String haplotypeName : haplotypeNames)
        {
            List<String> relevantEventIds = haplotypeNameToEventIds.get(haplotypeName);
            Map<String, Integer> eventIdToCount = relevantEventIds.stream().collect(Collectors.toMap(e -> e, e -> 2));

            Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(panel, eventIdToCount);
            HaplotypeCombination expectedBestHaplotypeCombination = new HaplotypeCombination(Map.of(haplotypeName, 2));

            assertEquals(1, geneToHaplotypeAnalysis.size());
            assertTrue(geneToHaplotypeAnalysis.containsKey(gene));
            assertEquals(expectedBestHaplotypeCombination, geneToHaplotypeAnalysis.get(gene).getBestHaplotypeCombination());
        }
    }

    private static void testAllSimpleHeterozygousHaplotypeCombinations(String gene, HaplotypePanel panel)
    {
        Map<String, List<String>> haplotypeNameToEventIds = panel.getNonDefaultHaplotypes(gene)
                .stream()
                .collect(Collectors.toMap(NonDefaultHaplotype::getName, PeachAlgoTest::getEventIds));
        String defaultHaplotypeName = panel.getDefaultHaplotype(gene).getName();
        haplotypeNameToEventIds.put(defaultHaplotypeName, Collections.emptyList());

        List<String> haplotypeNames = haplotypeNameToEventIds.keySet().stream().sorted().collect(Collectors.toList());
        for(int i = 0; i < haplotypeNames.size(); i++)
        {
            for(int j = i + 1; j < haplotypeNames.size(); j++)
            {
                String firstHaplotypeName = haplotypeNames.get(i);
                String secondHaplotypeName = haplotypeNames.get(j);

                Map<String, Integer> eventIdToCount = Stream.of(firstHaplotypeName, secondHaplotypeName)
                        .map(haplotypeNameToEventIds::get)
                        .flatMap(Collection::stream)
                        .collect(Collectors.groupingBy(Function.identity(), Collectors.collectingAndThen(Collectors.counting(), Long::intValue)));
                Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = determineGeneToHaplotypeAnalysis(panel, eventIdToCount);

                HaplotypeCombination expectedBestHaplotypeCombination =
                        new HaplotypeCombination(Map.of(firstHaplotypeName, 1, secondHaplotypeName, 1));

                assertEquals(1, geneToHaplotypeAnalysis.size());
                assertTrue(geneToHaplotypeAnalysis.containsKey(gene));
                assertEquals(expectedBestHaplotypeCombination, geneToHaplotypeAnalysis.get(gene).getBestHaplotypeCombination());
            }
        }
    }

    @NotNull
    private static List<String> getEventIds(final NonDefaultHaplotype h)
    {
        return h.events.stream().map(HaplotypeEvent::id).collect(Collectors.toList());
    }

    @NotNull
    private static Map<String, HaplotypeAnalysis> determineGeneToHaplotypeAnalysis(@NotNull HaplotypePanel haplotypePanel,
            @NotNull Map<String, Integer> eventIdToCount)
    {
        PeachAlgo algo = new PeachAlgo(haplotypePanel);
        return algo.getGeneToHaplotypeAnalysis(eventIdToCount);
    }

    private void assertEqualHaplotypeAnalysis(@NotNull HaplotypeAnalysis expected, @NotNull HaplotypeAnalysis actual)
    {
        assertEquals(expected.getDefaultHaplotypeName(), actual.getDefaultHaplotypeName());
        assertEquals(expected.getWildTypeHaplotypeName(), actual.getWildTypeHaplotypeName());

        assertEquals(expected.getEventIds(), actual.getEventIds());
        for(String eventId : expected.getEventIds())
        {
            assertEquals(String.format("Compare event counts of %s", eventId), expected.getEventCount(eventId), actual.getEventCount(eventId));
        }
        assertEqualHaplotypeCombinations(expected.getHaplotypeCombinations(), actual.getHaplotypeCombinations());
        assertEquals(expected.getQcStatus(), actual.getQcStatus());
        assertEquals(expected.getBestHaplotypeCombination(), actual.getBestHaplotypeCombination());
    }

    private void assertEqualHaplotypeCombinations(@NotNull List<HaplotypeCombination> expected, @NotNull List<HaplotypeCombination> actual)
    {
        assertEquals(expected.size(), actual.size());
        for(HaplotypeCombination expectedCombination : expected)
        {
            List<HaplotypeCombination> matchingActualCombinations =
                    actual.stream().filter(expectedCombination::equals).collect(Collectors.toList());
            assertEquals(String.format("Check exactly one combination matches %s", expectedCombination), 1, matchingActualCombinations.size());
        }
    }
}
