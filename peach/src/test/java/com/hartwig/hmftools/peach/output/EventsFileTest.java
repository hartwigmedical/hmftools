package com.hartwig.hmftools.peach.output;

import static junit.framework.TestCase.assertEquals;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

public class EventsFileTest
{
    private static final String EXPECTED_HEADER = "event\tcount";

    @Test
    public void testEmpty()
    {
        assertEquals(List.of(EXPECTED_HEADER), EventsFile.toLines(new HashMap<>()));
    }

    @Test
    public void testNonEmpty()
    {
        Map<String, Integer> eventIdToCount = new HashMap<>();
        eventIdToCount.put("EVENT_3", 1);
        eventIdToCount.put("EVENT_1", 2);
        eventIdToCount.put("EVENT_2", null);

        List<String> outputLines = EventsFile.toLines(eventIdToCount);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "EVENT_1\t2",
                "EVENT_2\tUNKNOWN",
                "EVENT_3\t1"
        );
        assertEquals(expectedLines, outputLines);
    }
}
