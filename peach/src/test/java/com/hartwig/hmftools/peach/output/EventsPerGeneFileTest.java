package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.PeachQCStatus;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.junit.Test;

public class EventsPerGeneFileTest
{
    private static final String EXPECTED_HEADER = "gene\tevent\tcount";

    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER), EventsPerGeneFile.toLines(new HashMap<>()));
    }

    @Test
    public void testNonEmpty()
    {
        HaplotypeAnalysis fake1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        HaplotypeAnalysis fake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_3", 1, "EVENT_1", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1))
                ),
                "*9",
                "*1",
                PeachQCStatus.FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND,
                null
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(
                        new HaplotypeCombination(Map.of("*4", 1, "*3", 1))
                ),
                "*9",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*4", 1, "*3", 1))
        );
        Map<String, Integer> fake4EventIdToCount = new HashMap<>();
        fake4EventIdToCount.put("EVENT_1", null);
        HaplotypeAnalysis fake4HaplotypeAnalysis = new HaplotypeAnalysis(
                fake4EventIdToCount,
                Collections.emptyList(),
                "*9",
                "*1",
                PeachQCStatus.FAIL_NO_COMBINATION_FOUND,
                null
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = Map.of(
                "FAKE4", fake4HaplotypeAnalysis,
                "FAKE3", fake3HaplotypeAnalysis,
                "FAKE2", fake2HaplotypeAnalysis,
                "FAKE1", fake1HaplotypeAnalysis
        );
        List<String> outputLines = EventsPerGeneFile.toLines(geneToHaplotypeAnalysis);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "FAKE2\tEVENT_1\t2",
                "FAKE2\tEVENT_3\t1",
                "FAKE3\tEVENT_1\t2",
                "FAKE4\tEVENT_1\tUNKNOWN"
        );
        assertEquals(expectedLines, outputLines);
    }
}
