package com.hartwig.hmftools.common.pileup;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PileupFileTest
{

    @Test
    public void testIndelSize() {
        final String text = ".$......+2AG.+3AGA.+11AAAAAAAAAAAAGGG";
        assertEquals(3, PileupFile.indelStringSize(9, text));
        assertEquals(4, PileupFile.indelStringSize(14, text));
        assertEquals(13, PileupFile.indelStringSize(20, text));
    }

    @Test
    public void testShortInsert() {
        final String line = "seq2\t156\tA\t11\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals("A", pileup.referenceBase());
        assertEquals(11, pileup.readCount());
        assertEquals(9, pileup.referenceCount());
        assertEquals(2, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(3, pileup.insertions());
        assertEquals(0, pileup.deletions());
    }

    @Test
    public void testLongInsert() {
        final String line = "seq2\t156\tA\t11\t.$......+2AG.+2AG.+12AGGGGGGGGGGGGG\t<975;:<<<<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(156, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(11, pileup.readCount());
        assertEquals(9, pileup.referenceCount());
        assertEquals(2, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(3, pileup.insertions());
        assertEquals(0, pileup.deletions());
    }

    @Test
    public void testShortDelete() {
        final String line = "seq3\t200\tA\t20\t,,,,,..,.-4CACC.-4CACC....,.,,.^~.\t==<<<<<<<<<<<::<;2<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(200, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(20, pileup.readCount());
        assertEquals(20, pileup.referenceCount());
        assertEquals(0, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(0, pileup.insertions());
        assertEquals(2, pileup.deletions());
    }

    @Test
    public void testLongDelete() {
        final String line = "seq3\t200\tA\t20\t,,,,,..,.-13CACCCCCCCCCCCC-4CACC....,.,,.^~.\t==<<<<<<<<<<<::<;2<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(200, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(20, pileup.readCount());
        assertEquals(19, pileup.referenceCount());
        assertEquals(0, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(1, pileup.cMismatchCount());
        assertEquals(0, pileup.insertions());
        assertEquals(2, pileup.deletions());
    }
}
