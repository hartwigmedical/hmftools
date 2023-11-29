package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.junit.Test;

public class BestHaplotypeCombinationsFileTest
{
    private static final String EXPECTED_HEADER = "gene\tallele\tcount";
    
    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER),  BestHaplotypeCombinationsFile.toLines(new HashMap<>()));
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
                Map.of("EVENT_1", 1, "EVENT_3", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1))
                ),
                "*9",
                "*1"
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(
                        new HaplotypeCombination(Map.of("*4", 1, "*3", 1))
                ),
                "*9",
                "*1"
        );
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = Map.of(
                "FAKE3", fake3HaplotypeAnalysis,
                "FAKE2", fake2HaplotypeAnalysis,
                "FAKE1", fake1HaplotypeAnalysis
        );
        List<String> outputLines = BestHaplotypeCombinationsFile.toLines(geneToHaplotypeAnalysis);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "FAKE1\t*1\t2",
                "FAKE2\tUNRESOLVED\t2",
                "FAKE3\t*3\t1",
                "FAKE3\t*4\t1"
        );
        assertEquals(expectedLines, outputLines);
    }
}
