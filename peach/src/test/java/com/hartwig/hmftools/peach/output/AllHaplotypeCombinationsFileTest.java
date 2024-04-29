package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.PeachQCStatus;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

public class AllHaplotypeCombinationsFileTest
{
    private static final String EXPECTED_HEADER = "gene\tcombination\tcount\tnonWildTypeCount";

    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER), AllHaplotypeCombinationsFile.toLines(new HashMap<>()));
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
                Map.of("EVENT_1", 1, "EVENT_3", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*1", 2)),
                        new HaplotypeCombination(Map.of("*3", 1, "*1", 1))
                ),
                "*9",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*3", 1, "*1", 1))
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(new HaplotypeCombination(Map.of("*9", 2))),
                "*9",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*9", 2))
        );
        HaplotypeAnalysis fake4HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(),
                "*9",
                "*1",
                PeachQCStatus.FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND,
                null
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = Map.of(
                "FAKE3", fake3HaplotypeAnalysis,
                "FAKE2", fake2HaplotypeAnalysis,
                "FAKE1", fake1HaplotypeAnalysis
        );
        List<String> outputLines = AllHaplotypeCombinationsFile.toLines(geneToHaplotypeAnalysis);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "FAKE1\t(*1,2)\t2\t0",
                "FAKE2\t(*1,1);(*3,1)\t2\t1",
                "FAKE2\t(*9,1);(2373C>T,1)\t2\t2",
                "FAKE2\t(*1,2);(*2,1)\t3\t1",
                "FAKE3\t(*9,2)\t2\t2"
        );
        assertEquals(expectedLines, outputLines);
    }
}
