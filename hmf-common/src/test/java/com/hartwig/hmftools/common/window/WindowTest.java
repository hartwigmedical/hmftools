package com.hartwig.hmftools.common.window;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.window.Window;

import org.junit.Test;

public class WindowTest {

    @Test
    public void testWindowStart() {
        final Window victim = new Window(1000);

        assertEquals(1, victim.start(1));
        assertEquals(1, victim.start(1000));

        assertEquals(1001, victim.start(1001));
        assertEquals(1001, victim.start(2000));

        assertEquals(2001, victim.start(2001));
        assertEquals(2001, victim.start(3000));

        assertEquals(249250001, victim.start(249250621));
    }

    @Test
    public void testWindowEnd() {
        final Window victim = new Window(1000);

        assertEquals(1000, victim.end(1));
        assertEquals(1000, victim.end(1000));

        assertEquals(2000, victim.end(1001));
        assertEquals(2000, victim.end(2000));

        assertEquals(3000, victim.end(2001));
        assertEquals(3000, victim.end(3000));

        assertEquals(249251000, victim.end(249250621));
    }
}
