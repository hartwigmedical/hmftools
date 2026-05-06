package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class WindowTest
{
    @Test
    public void testWindowStart()
    {
        Window window = new Window(1000);

        assertEquals(1, window.start(1));
        assertEquals(1, window.start(1000));

        assertEquals(1001, window.start(1001));
        assertEquals(1001, window.start(2000));

        assertEquals(2001, window.start(2001));
        assertEquals(2001, window.start(3000));

        assertEquals(249250001, window.start(249250621));
    }

    @Test
    public void testWindowEnd()
    {
        Window window = new Window(1000);

        assertEquals(1000, window.end(1));
        assertEquals(1000, window.end(1000));

        assertEquals(2000, window.end(1001));
        assertEquals(2000, window.end(2000));

        assertEquals(3000, window.end(2001));
        assertEquals(3000, window.end(3000));

        assertEquals(249251000, window.end(249250621));
    }
}
