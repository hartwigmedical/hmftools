package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ReadContextImprovedTest {

    private ReadContextImproved victim = new ReadContextImproved(Strings.EMPTY, "", 1000, 5, 4, 6, 3, "GATCTCCTCA".getBytes());

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

        assertEquals(ReadContextMatch.FULL, victim.centreMatch(5, "GATCTCCTCA".getBytes()));
        assertEquals(ReadContextMatch.FULL, victim.centreMatch(1, "TCC".getBytes()));

        assertEquals(ReadContextMatch.NONE, victim.centreMatch(1, "CCC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(1, "TTC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(1, "TCT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(1, "TC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "CC".getBytes()));
    }

    @Test
    public void testPloyAJitterCentreMatch() {

        ReadContextImproved victim = new ReadContextImproved(Strings.EMPTY, "A", 1000, 2, 2, 11, 2, "GGGAAAAAAAATTT".getBytes());
        assertEquals(ReadContextMatch.FULL, victim.centreMatch(0, "GAAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAACAAAT".getBytes()));

        assertEquals(ReadContextMatch.JITTER_ADDED, victim.centreMatch(0, "GAAAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAATAAAAT".getBytes()));

        assertEquals(ReadContextMatch.JITTER_REMOVED, victim.centreMatch(0, "GAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAGAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAA".getBytes()));


        victim = new ReadContextImproved(Strings.EMPTY, "", 1000, 2, 2, 11, 2, "GGGAAAAAAAATTT".getBytes());
        assertEquals(ReadContextMatch.FULL, victim.centreMatch(0, "GAAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAACAAAT".getBytes()));

        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAATAAAAT".getBytes()));

        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAAAAC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAAAGAAT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GAAA".getBytes()));

    }

    @Test
    public void testRepeatJitterCentreMatch() {

         ReadContextImproved victim = new ReadContextImproved(Strings.EMPTY, "AC", 1000, 2, 2, 11, 2, "GGGACACACACTTT".getBytes());
        assertEquals(ReadContextMatch.FULL, victim.centreMatch(0, "GACACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACCCT".getBytes()));

        assertEquals(ReadContextMatch.JITTER_ADDED, victim.centreMatch( 0, "GACACACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACACCCT".getBytes()));

        assertEquals(ReadContextMatch.JITTER_REMOVED, victim.centreMatch( 0, "GACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACTCT".getBytes()));

        victim = new ReadContextImproved(Strings.EMPTY, "", 1000, 2, 2, 11, 2, "GGGACACACACTTT".getBytes());
        assertEquals(ReadContextMatch.FULL, victim.centreMatch(0, "GACACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACCCT".getBytes()));

        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch(0, "GACACACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACACCCT".getBytes()));

        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACACC".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.centreMatch( 0, "GACACTCT".getBytes()));
    }

    @Test
    public void testPartialMatchMustHaveAtLeastOneFullSide() {

        ReadContextImproved victim = new ReadContextImproved(Strings.EMPTY, "", 1000, 2, 2, 2, 2, "GGTAA".getBytes());
        assertEquals(ReadContextMatch.FULL, victim.matchAtPosition(2, "GGTAA".getBytes()));

        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(2, "GGTA".getBytes()));
        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(2, "GGT".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.matchAtPosition(1, "GT".getBytes()));

        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(1, "GTAA".getBytes()));
        assertEquals(ReadContextMatch.PARTIAL, victim.matchAtPosition(0, "TAA".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.matchAtPosition(0, "TA".getBytes()));
        assertEquals(ReadContextMatch.NONE, victim.matchAtPosition(0, "T".getBytes()));
    }





}
