package com.hartwig.hmftools.cobalt.count;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChromosomeReadCountTest {

    @Test
    public void testWindowPosition() {
        final int windowSize = 1000;

        assertEquals(1, ChromosomeReadCount.windowPosition(windowSize, 1));
        assertEquals(1, ChromosomeReadCount.windowPosition(windowSize, 1000));

        assertEquals(1001, ChromosomeReadCount.windowPosition(windowSize, 1001));
        assertEquals(1001, ChromosomeReadCount.windowPosition(windowSize, 2000));

        assertEquals(2001, ChromosomeReadCount.windowPosition(windowSize, 2001));
        assertEquals(2001, ChromosomeReadCount.windowPosition(windowSize, 3000));

        assertEquals(249250001, ChromosomeReadCount.windowPosition(windowSize, 249250621));
    }

}
