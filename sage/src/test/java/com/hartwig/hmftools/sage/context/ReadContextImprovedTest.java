package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ReadContextImprovedTest {

    private ReadContextImproved victim = new ReadContextImproved(1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

    @Test
    public void testRightFlankMatchingBases() {
        assertEquals(-1, victim.rightFlankMatchingBases(3, "TCTCCTCG".getBytes()));

        assertEquals(3, victim.rightFlankMatchingBases(3, "TCTCCTCAG".getBytes()));
        assertEquals(3, victim.rightFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(2, victim.rightFlankMatchingBases(3, "TCTCCTC".getBytes()));
        assertEquals(1, victim.rightFlankMatchingBases(3, "TCTCCT".getBytes()));
        assertEquals(0, victim.rightFlankMatchingBases(3, "TCTCC".getBytes()));
    }

    @Test
    public void testLeftFlankMatchingBases() {

        assertEquals(-1, victim.leftFlankMatchingBases(5, "TTCTCCTCA".getBytes()));

        assertEquals(3, victim.leftFlankMatchingBases(5, "GATCTCCTCA".getBytes()));
        assertEquals(3, victim.leftFlankMatchingBases(4, "ATCTCCTCA".getBytes()));
        assertEquals(2, victim.leftFlankMatchingBases(3, "TCTCCTCA".getBytes()));
        assertEquals(1, victim.leftFlankMatchingBases(2, "CTCCTCA".getBytes()));
        assertEquals(0, victim.leftFlankMatchingBases(1, "TCCTCA".getBytes()));
    }

    @Test
    public void testCentreMatch() {

        assertTrue(victim.centerMatch(5, "GATCTCCTCA".getBytes()));
        assertTrue(victim.centerMatch(1, "TCC".getBytes()));

        assertFalse(victim.centerMatch(1, "CCC".getBytes()));
        assertFalse(victim.centerMatch(1, "TTC".getBytes()));
        assertFalse(victim.centerMatch(1, "TCT".getBytes()));
        assertFalse(victim.centerMatch(1, "TC".getBytes()));
        assertFalse(victim.centerMatch(0, "CC".getBytes()));
    }

    @Test
    public void testPloyAJitterCentreMatch() {

        final ReadContextImproved victim = new ReadContextImproved(1000, 2, 2, 11, 2, "GGGAAAAAAAATTT".getBytes());
        assertTrue(victim.centerMatch(0, 0, "GAAAAAAAAT".getBytes()));
        assertFalse(victim.centerMatch(0, 0, "GAAAAAAAAC".getBytes()));
        assertFalse(victim.centerMatch(0, 0, "GAAAACAAAT".getBytes()));

        assertTrue(victim.centerMatch(1, 0, "GAAAAAAAAAT".getBytes()));
        assertFalse(victim.centerMatch(1, 0, "GAAAAAAAAAC".getBytes()));
        assertFalse(victim.centerMatch(1, 0, "GAAAATAAAAT".getBytes()));

        assertTrue(victim.centerMatch(-1, 0, "GAAAAAAAT".getBytes()));
        assertFalse(victim.centerMatch(-1, 0, "GAAAAAAAC".getBytes()));
        assertFalse(victim.centerMatch(-1, 0, "GAAAAGAAT".getBytes()));
    }

    @Test
    public void testRepeatJitterCentreMatch() {

        final ReadContextImproved victim = new ReadContextImproved(1000, 2, 2, 11, 2, "GGGACACACACTTT".getBytes());
        assertTrue(victim.centerMatch(0, 0, "GACACACACT".getBytes()));
        assertFalse(victim.centerMatch(0, 0, "GACACACACC".getBytes()));
        assertFalse(victim.centerMatch(0, 0, "GACACACCCT".getBytes()));

        assertTrue(victim.centerMatch(2, 0, "GACACACACACT".getBytes()));
        assertFalse(victim.centerMatch(2, 0, "GACACACACACC".getBytes()));
        assertFalse(victim.centerMatch(2, 0, "GACACACACCCT".getBytes()));

        assertTrue(victim.centerMatch(-2, 0, "GACACACT".getBytes()));
        assertFalse(victim.centerMatch(-2, 0, "GACACACC".getBytes()));
        assertFalse(victim.centerMatch(-2, 0, "GACACTCT".getBytes()));
    }

}
