package com.hartwig.hmftools.serve.extraction.events;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class EventInterpretationFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String eventInterpretationTsv =
                EventInterpretationFile.eventInterpretationTsv(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR);
        List<EventInterpretation> eventInterpretations = EventInterpretationFile.read(eventInterpretationTsv);

        assertEquals(2, eventInterpretations.size());

        List<String> lines = EventInterpretationFile.toLines(eventInterpretations);
        List<EventInterpretation> interpretationEvents = EventInterpretationFile.fromLines(lines);
        List<String> regeneratedLines = EventInterpretationFile.toLines(interpretationEvents);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }

}