package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReadContextTest {

    @Test
    public void testExactLengthOnBothSide() {
        final String expected = "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC";
        final ReadContext victim = new ReadContext(15, 15, expected.getBytes());
        assertTrue(victim.isComplete());
        assertEquals(expected, victim.toString());
    }

    @Test
    public void testAdditionalLengthOnBothSide() {
        final String expected = "ATG" + "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC" + "ATG";
        final ReadContext victim = new ReadContext(15, 18, expected.getBytes());
        assertTrue(victim.isComplete());
        assertEquals("AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC", victim.toString());
    }

    @Test
    public void testShortLeft() {
        final String bytes = "AATCCCCCCCCCCCCCCC";
        final ReadContext victim = new ReadContext(15, 2, bytes.getBytes());
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testShortRight() {
        final String bytes = "AAAAAAAAAAAAAAATCCC";
        final ReadContext victim = new ReadContext(15, 15, bytes.getBytes());
        assertFalse(victim.isComplete());
        assertEquals(bytes, victim.toString());
    }

    @Test
    public void testIncompleteStillMatches() {
        final ReadContext full = new ReadContext(15, 15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes());
        final ReadContext shortLeft = new ReadContext(15, 2, "AATCCCCCCCCCCCCCCC".getBytes());
        final ReadContext shortRight = new ReadContext(15, 15, "AAAAAAAAAAAAAAATCCC".getBytes());

        assertEquals(ReadContextMatch.FULL, full.match(full));

        assertEquals(ReadContextMatch.PARTIAL, full.match(shortLeft));
        assertEquals(ReadContextMatch.PARTIAL, full.match(shortRight));

        assertEquals(ReadContextMatch.PARTIAL, shortLeft.match(full));
        assertEquals(ReadContextMatch.PARTIAL, shortRight.match(full));

        assertEquals(ReadContextMatch.NONE, shortLeft.match(shortRight));
        assertEquals(ReadContextMatch.NONE, shortRight.match(shortLeft));
    }

    @Test
    public void testMatch() {
        final ReadContext full = new ReadContext(15, 15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes());
        final ReadContext diffAtAlt = new ReadContext(15, 15, "AAAAAAAAAAAAAAAGCCCCCCCCCCCCCCC".getBytes());
        final ReadContext diffLeft = new ReadContext(15, 15, "TAAAAAAAAAAAAAATCCCCCCCCCCCCCCC".getBytes());
        final ReadContext diffRight = new ReadContext(15, 15, "AAAAAAAAAAAAAAATCCCCCCCCCCCCCCG".getBytes());
        final ReadContext diffBoth = new ReadContext(15, 15, "TAAAAAAAAAAAAAATCCCCCCCCCCCCCCG".getBytes());

        assertEquals(ReadContextMatch.NONE, full.match(diffAtAlt));
        assertEquals(ReadContextMatch.NONE, full.match(diffLeft));
        assertEquals(ReadContextMatch.NONE, full.match(diffRight));
        assertEquals(ReadContextMatch.NONE, full.match(diffBoth));
    }

}
