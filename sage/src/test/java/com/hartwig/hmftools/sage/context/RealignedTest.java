package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.context.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.context.RealignedType.NONE;
import static com.hartwig.hmftools.sage.context.RealignedType.SHORTENED;

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
        int endIndex = startIndex + 55;
        Realigned victim = new Realigned();

        assertRealigned(EXACT, 0, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequence.getBytes()));
        assertRealigned(EXACT, 0, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequenceMinusTA.getBytes()));
        assertRealigned(EXACT, 0, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequencePlusTA.getBytes()));
        assertRealigned(EXACT, 0, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequenceMinusGT.getBytes()));
        assertRealigned(EXACT, 0, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequencePlusGT.getBytes()));

        assertRealigned(SHORTENED, 9, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequenceMinusTA.getBytes()));
        assertRealigned(LENGTHENED, 11, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequencePlusTA.getBytes()));

        assertRealigned(SHORTENED, 6, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequenceMinusGT.getBytes()));
        assertRealigned(LENGTHENED, 8, victim.realigned(startIndex, endIndex, sequence.getBytes(), sequencePlusGT.getBytes()));

        assertRealigned(SHORTENED, 10, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequence.getBytes()));
        assertRealigned(LENGTHENED, 10, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequence.getBytes()));

        assertRealigned(SHORTENED, 7, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequence.getBytes()));
        assertRealigned(LENGTHENED, 7, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequence.getBytes()));

        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequencePlusTA.getBytes(), sequenceMinusTA.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequenceMinusTA.getBytes(), sequencePlusTA.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequenceMinusGT.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequencePlusGT.getBytes()));

        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequencePlusTA.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequenceMinusGT.getBytes(), sequenceMinusTA.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequencePlusTA.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(startIndex, endIndex, sequencePlusGT.getBytes(), sequenceMinusTA.getBytes()));
    }

    @Test
    public void testPolyA() {
        String sequence = "GAAAAAT";
        String sequencePlusA = "GAAAAAAT";
        String sequenceMinusA = "GAAAAT";

        Realigned victim = new Realigned();
        assertRealigned(EXACT, 0, victim.realigned(0, 6, sequence.getBytes(), sequence.getBytes()));
        assertRealigned(LENGTHENED, 6, victim.realigned(0, 6, sequence.getBytes(), sequencePlusA.getBytes()));
        assertRealigned(SHORTENED, 4, victim.realigned(0, 6, sequence.getBytes(), sequenceMinusA.getBytes()));

        assertRealigned(LENGTHENED, 5, victim.realigned(0, 5, sequenceMinusA.getBytes(), sequence.getBytes()));
        assertRealigned(SHORTENED, 5, victim.realigned(0, 7, sequencePlusA.getBytes(), sequence.getBytes()));

    }

    @Test
    public void testMinRepeatCount() {
        String sequence3A = "GAAAT";
        String sequence4A = "GAAAAT";
        String sequence5A = "GAAAAAT";

        Realigned victim = new Realigned();
        assertRealigned(EXACT, 0, victim.realigned(0, 4, sequence3A.getBytes(), sequence3A.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(0, 4, sequence3A.getBytes(), sequence4A.getBytes()));
        assertRealigned(NONE, 0, victim.realigned(0, 5, sequence4A.getBytes(), sequence3A.getBytes()));
        assertRealigned(LENGTHENED, 5, victim.realigned(0, 5, sequence4A.getBytes(), sequence5A.getBytes()));
        assertRealigned(SHORTENED, 4, victim.realigned(0, 6, sequence5A.getBytes(), sequence4A.getBytes()));
    }

    @Test
    public void testReadRepeatTakesReadContextIndexBeforeZero() {
        String sequence = "ATACTAAAAAAAAAAAAAAAAAAAA";
        String read = "GATACGATACGATACGATACGATACCAAAAAAAAAA";

        Realigned victim = new Realigned();
        assertRealigned(NONE, 0, victim.realigned(0, 4, sequence.getBytes(), read.getBytes()));
    }

    private static void assertRealigned(RealignedType expectedType, int expectedCount, RealignedContext context) {
        assertEquals(expectedCount, context.repeatCount());
        assertEquals(expectedType, context.type());
    }

}
