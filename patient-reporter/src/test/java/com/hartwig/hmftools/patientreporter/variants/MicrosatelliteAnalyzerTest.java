package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.repeat.ImmutableRepeatContext;
import com.hartwig.hmftools.common.purple.repeat.RepeatContext;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MicrosatelliteAnalyzerTest {

    @Test
    public void testShortRepeatContextRelevance() {
        assertRepeatContextRelevance(false, 4, "A");
        assertRepeatContextRelevance(true, 5, "A");
        assertRepeatContextRelevance(true, 100, "A");
    }

    @Test
    public void testLongRepeatContextRelevance() {
        assertRepeatContextRelevance(false, 3, "AT");
        assertRepeatContextRelevance(true, 4, "AT");
        assertRepeatContextRelevance(true, 100, "AT");

        assertRepeatContextRelevance(false, 3, "ATG");
        assertRepeatContextRelevance(true, 4, "ATG");
        assertRepeatContextRelevance(true, 100, "ATG");

        assertRepeatContextRelevance(false, 3, "ATGA");
        assertRepeatContextRelevance(true, 4, "ATGA");
        assertRepeatContextRelevance(true, 100, "ATGA");

        assertRepeatContextRelevance(false, 3, "ATGAA");
        assertRepeatContextRelevance(false, 4, "ATGAA");
    }

    private static void assertRepeatContextRelevance(boolean expectedResult, int count, @NotNull String sequence) {
        final RepeatContext context = ImmutableRepeatContext.builder().count(count).sequence(sequence).build();
        assertEquals(expectedResult, MicrosatelliteAnalyzer.repeatContextIsRelevant(context));
    }
}
