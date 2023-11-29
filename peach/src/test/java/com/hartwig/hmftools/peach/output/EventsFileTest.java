package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.junit.Test;

public class EventsFileTest
{
    private static final String EXPECTED_HEADER = "event\tcount";
    
    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER),  EventsFile.toLines(new HashMap<>()));
    }

    @Test
    public void testNonEmpty()
    {
        Map<String, Integer> eventIdToCount = Map.of("EVENT_3", 1, "EVENT_1", 2);
        List<String> outputLines = EventsFile.toLines(eventIdToCount);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "EVENT_1\t2",
                "EVENT_3\t1"
        );
        assertEquals(expectedLines, outputLines);
    }
}
