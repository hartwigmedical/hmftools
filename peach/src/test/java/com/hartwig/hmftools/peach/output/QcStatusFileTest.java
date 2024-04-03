package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.junit.Test;

public class QcStatusFileTest
{
    private static final String EXPECTED_HEADER = "gene\tstatus";

    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER), QcStatusFile.toLines(new HashMap<>()));
    }

    @Test
    public void testNonEmpty()
    {
        HaplotypeAnalysis fake1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1"
        );
        HaplotypeAnalysis fake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_3", 1, "EVENT_1", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1))
                ),
                "*9",
                "*1"
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                new ArrayList<>(),
                "*9",
                "*1"
        );
        HaplotypeAnalysis fake4HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_5", 1, "EVENT_6", 2),
                List.of(new HaplotypeCombination(Map.of("*3", 2, "*2", 1))),
                "*9",
                "*1"
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = Map.of(
                "FAKE3", fake3HaplotypeAnalysis,
                "FAKE2", fake2HaplotypeAnalysis,
                "FAKE1", fake1HaplotypeAnalysis,
                "FAKE4", fake4HaplotypeAnalysis
        );
        List<String> outputLines = QcStatusFile.toLines(geneToHaplotypeAnalysis);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "FAKE1\tPASS",
                "FAKE2\tFAIL_NO_UNIQUE_BEST_COMBINATION_FOUND",
                "FAKE3\tFAIL_NO_COMBINATION_FOUND",
                "FAKE4\tWARN_TOO_MANY_ALLELES_FOUND"
        );
        assertEquals(expectedLines, outputLines);
    }
}