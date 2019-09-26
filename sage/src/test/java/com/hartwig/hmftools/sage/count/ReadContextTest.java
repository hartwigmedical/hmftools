package com.hartwig.hmftools.sage.count;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReadContextTest {

    @Test
    public void testSufficientLengthOnBothSide() {
        final String bytes = "AAAAAAAAAATCCCCCCCCCC";
        final ReadContext victim = new ReadContext(10, bytes.getBytes());
        final ReadContext expected = new ReadContext("AAAAAAAAAA", "T", "CCCCCCCCCC");
        assertEquals(expected.toString(), victim.toString());
    }

    @Test
    public void testShortLeft() {
        final String bytes = "AATCCCCCCCCCC";
        final ReadContext victim = new ReadContext(3, bytes.getBytes());
        final ReadContext expected = new ReadContext("AA", "T", "CCCCCCCCCC");
        assertEquals(expected.toString(), victim.toString());
    }

    @Test
    public void testShortRight() {
        final String bytes = "AAAAAAAAAATCCC";
        final ReadContext victim = new ReadContext(10, bytes.getBytes());
        final ReadContext expected = new ReadContext("AAAAAAAAAA", "T", "CCC");
        assertEquals(expected.toString(), victim.toString());
    }

    @Test
    public void testIncompleteStillMatches() {
        final ReadContext full = new ReadContext("AAAAAAAAAA", "T", "CCCCCCCCCC");
        final ReadContext shortLeft = new ReadContext("AA", "T", "CCCCCCCCCC");
        final ReadContext shortRight = new ReadContext("AAAAAAAAAA", "T", "CCC");

        assertTrue(full.match(shortLeft));
        assertTrue(full.match(shortRight));

        assertTrue(shortLeft.match(full));
        assertTrue(shortRight.match(full));

        assertFalse(shortLeft.match(shortRight));
        assertFalse(shortRight.match(shortLeft));
    }

}
