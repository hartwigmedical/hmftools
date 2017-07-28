package com.hartwig.hmftools.cobalt;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChromosomeCountTest {

    @Test
    public void testWindow() {
        final ChromosomeCount count = new ChromosomeCount("chr1", "1", 249250621, 1000);
        assertEquals(1, count.windowPosition(1));
        assertEquals(1, count.windowPosition(1000));

        assertEquals(1001, count.windowPosition(1001));
        assertEquals(1001, count.windowPosition(2000));

        assertEquals(2001, count.windowPosition(2001));
        assertEquals(2001, count.windowPosition(3000));
    }

    @Test
    public void testFinalWindow() {
        final ChromosomeCount count = new ChromosomeCount("chr1", "1", 249250621, 1000);
        assertEquals(249250001, count.lastWindowPosition());
    }

}
