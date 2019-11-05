package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.context.MatchType;

import org.apache.logging.log4j.util.Strings;
import org.junit.Assert;
import org.junit.Test;

public class ReadContextTest {

    private ReadContext victim = new ReadContext("", 1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

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

        assertEquals(true, victim.centreMatch(5, "GATCTCCTCA".getBytes()));
        assertEquals(true, victim.centreMatch(1, "TCC".getBytes()));

        assertEquals(false, victim.centreMatch(1, "CCC".getBytes()));
        assertEquals(false, victim.centreMatch(1, "TTC".getBytes()));
        assertEquals(false, victim.centreMatch(1, "TCT".getBytes()));
        assertEquals(false, victim.centreMatch(1, "TC".getBytes()));
        assertEquals(false, victim.centreMatch(0, "CC".getBytes()));
    }


    @Test
    public void testPartialMatchMustHaveAtLeastOneFullSide() {

        ReadContext victim = new ReadContext("", 1000, 2, 2, 2, 2, "GGTAA".getBytes());
        Assert.assertEquals(MatchType.FULL, victim.matchAtPosition(2, "GGTAA".getBytes()));

        assertEquals(MatchType.PARTIAL, victim.matchAtPosition(2, "GGTA".getBytes()));
        assertEquals(MatchType.PARTIAL, victim.matchAtPosition(2, "GGT".getBytes()));
        assertEquals(MatchType.NONE, victim.matchAtPosition(1, "GT".getBytes()));

        assertEquals(MatchType.PARTIAL, victim.matchAtPosition(1, "GTAA".getBytes()));
        assertEquals(MatchType.PARTIAL, victim.matchAtPosition(0, "TAA".getBytes()));
        assertEquals(MatchType.NONE, victim.matchAtPosition(0, "TA".getBytes()));
        assertEquals(MatchType.NONE, victim.matchAtPosition(0, "T".getBytes()));
    }

    @Test
    public void testNegativeReadIndex() {
        ReadContext victim = new ReadContext("", 1000, 2, 2, 2, 2, "GGTAA".getBytes());
        assertEquals(MatchType.FULL, victim.matchAtPosition(2, "GGTAA".getBytes()));
        assertEquals(MatchType.NONE, victim.matchAtPosition(-1, "GGTAA".getBytes()));
    }

    @Test
    public void testPhasedMNV() {
        ReadContext victim1 = new ReadContext(Strings.EMPTY, 1000, 4, 4, 4, 4, "GATCTTGATC".getBytes());
        ReadContext victim2 = new ReadContext(Strings.EMPTY, 1001, 4, 4, 5, 4, "ATCTTGATCT".getBytes());

        assertTrue(victim1.phased(victim2));
        assertTrue(victim2.phased(victim1));
    }

}
