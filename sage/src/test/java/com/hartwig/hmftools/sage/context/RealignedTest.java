package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RealignedTest {

    @Test
    public void testDiNucleotideRealigned() {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        String sequenceMinusTA = "CAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATCACATTTTTATTATTC";
        String sequencePlusTA = "TAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATATCACATTTTTATTATTT";
        String sequenceMinusGT = "GAGAGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        String sequencePlusGT = "GAGAGTGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";

        int startIndex = 3;
        int endIndex = startIndex + 50;
        Realigned victim = new Realigned();

        assertEquals(RealignedType.EXACT, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequence.getBytes()));
        assertEquals(RealignedType.EXACT, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequenceMinusTA.getBytes()));
        assertEquals(RealignedType.EXACT, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequencePlusTA.getBytes()));
        assertEquals(RealignedType.EXACT, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequenceMinusGT.getBytes()));
        assertEquals(RealignedType.EXACT, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequencePlusGT.getBytes()));

        assertEquals(RealignedType.SHORTENED, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequenceMinusTA.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequencePlusTA.getBytes()));

        assertEquals(RealignedType.SHORTENED, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequenceMinusGT.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequencePlusGT.getBytes()));

        assertEquals(RealignedType.SHORTENED, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequence.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequence.getBytes()));

        assertEquals(RealignedType.SHORTENED, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequence.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequence.getBytes()));

        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequenceMinusTA.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequencePlusTA.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequenceMinusGT.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequencePlusGT.getBytes()));

        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequencePlusTA.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequenceMinusTA.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequencePlusTA.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequenceMinusTA.getBytes()));

    }

    @Test
    public void testPolyA() {
        String sequence = "GAAAAAT";
        String sequencePlusA = "GAAAAAAT";
        String sequenceMinusA = "GAAAAT";

        Realigned victim = new Realigned();
        assertEquals(RealignedType.EXACT, victim.realigned(0, 6, sequence.getBytes(), sequence.getBytes()));
        assertEquals(RealignedType.SHORTENED, victim.realigned(0, 6, sequence.getBytes(), sequenceMinusA.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(0, 6, sequence.getBytes(), sequencePlusA.getBytes()));
    }

    @Test
    public void testMinRepeatCount() {
        String sequence3A = "GAAAT";
        String sequence4A = "GAAAAT";
        String sequence5A = "GAAAAAT";

        Realigned victim = new Realigned();
        assertEquals(RealignedType.EXACT, victim.realigned(0, 4, sequence3A.getBytes(), sequence3A.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(0, 4, sequence3A.getBytes(), sequence4A.getBytes()));
        assertEquals(RealignedType.NONE, victim.realigned(0, 5, sequence4A.getBytes(), sequence3A.getBytes()));
        assertEquals(RealignedType.LENGTHENED, victim.realigned(0, 5, sequence4A.getBytes(), sequence5A.getBytes()));
        assertEquals(RealignedType.SHORTENED, victim.realigned(0, 6, sequence5A.getBytes(), sequence4A.getBytes()));

    }

    @Test
    public void testReadRepeatTakesReadContextIndexBeforeZero() {
        String sequence = "ATACTAAAAAAAAAAAAAAAAAAAA";
        String read = "GATACGATACGATACGATACGATACCAAAAAAAAAA";

        Realigned victim = new Realigned();
        assertEquals(RealignedType.NONE, victim.realigned(0, 4, sequence.getBytes(), read.getBytes()));
    }

}
