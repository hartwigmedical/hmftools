package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;
import static com.hartwig.hmftools.peach.panel.TestHaplotypePanelFactory.createDefaultTestHaplotypePanel;

import static junit.framework.TestCase.assertEquals;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HaplotypeCallerTest
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
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1"
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
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1"
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
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CAT", 1,
                        "VAR_chr2_234668879_CAT_C", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*28", 1, "*36", 1))),
                "*1",
                "*1"
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
                Map.of(
                        "VAR_chr1_98348885_G_A", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*9A", 1, "*1", 1))),
                "*9A",
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234669144_G_A", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*1", 1, "*6", 1))),
                "*1",
                "*1"
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
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CATAT", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*1", 1, "*37", 1))),
                "*1",
                "*1"
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
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr2_234668879_C_CAT", 2
                ),
                List.of(new HaplotypeCombination(Map.of("*28", 2))),
                "*1",
                "*1"
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
                "*1"
        );
        HaplotypeAnalysis expectedFake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_C_CAT", 1
                ),
                List.of(new HaplotypeCombination(Map.of("*2", 1, "*1", 1))),
                "*1",
                "*1"
        );
        HaplotypeAnalysis expectedFake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of(
                        "VAR_chr3_1234567_CAT_C", 1,
                        "VAR_chr3_1234567_C_CAT", 1
                ),
                new ArrayList<>(),
                "*1",
                "*1"
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
    public void testAmbiguousCallingPreferFewerHaplotypes()
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
                "*1"
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
        assertTrue(expectedHaplotypeAnalysis.hasBestHaplotypeCombination());
        assertEquals(new HaplotypeCombination(Map.of("*2", 1, "*3", 1)), expectedHaplotypeAnalysis.getBestHaplotypeCombination());
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
                "*1"
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
        assertFalse(expectedHaplotypeAnalysis.hasBestHaplotypeCombination());
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
                "*1"
        );
        assertEquals(1, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("FAKE1"));

        assertEqualHaplotypeAnalysis(expectedHaplotypeAnalysis, geneToHaplotypeAnalysis.get("FAKE1"));
        assertTrue(expectedHaplotypeAnalysis.hasBestHaplotypeCombination());
        assertEquals(new HaplotypeCombination(Map.of("*2", 1, "*1", 1)), expectedHaplotypeAnalysis.getBestHaplotypeCombination());
    }

    @NotNull
    private static Map<String, HaplotypeAnalysis> determineGeneToHaplotypeAnalysis(@NotNull HaplotypePanel haplotypePanel,
            @NotNull Map<String, Integer> eventIdToCount)
    {
        HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
        return caller.getGeneToHaplotypeAnalysis(eventIdToCount);
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
