package com.hartwig.hmftools.common.pileup;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PileupFileTest {

    @Test
    public void testIndelSize() {
        final String text = ".$......+2AG.+3AGA.+11AAAAAAAAAAAAGGG";
        assertEquals(2, PileupFile.indelSize(9, text));
        assertEquals(3, PileupFile.indelSize(14, text));
        assertEquals(11, PileupFile.indelSize(20, text));
    }

    @Test
    public void testAltIndel() {
        final String line =
                "6\t106556396\tC\t73\t,.,,,,,.,..,.,,,,+1t.+1T..,....,+1t,+1t.,+1t.,+1t,+1t...,+3ttt,-1t..,+2tt.......+1T.+1T,+1t,+1t,+2tt,+1tt+1t..+1T,-2tt.+1T.+1T.+2TT,,+1t.t+1t...+2TT.+1T,-3ttt...-1T.+1T.";
        final Pileup pileup = PileupFile.fromString(line);

        assertEquals("C", pileup.referenceBase());
        assertEquals(73, pileup.readCount());
        assertEquals(44, pileup.referenceCount());
        assertEquals(0, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(2, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(25, pileup.insertions());
        assertEquals(4, pileup.deletions());
        assertEquals(29, pileup.indels());
        assertEquals(1, pileup.inframeInsertions());
        assertEquals(1, pileup.inframeDeletions());
        assertEquals(2, pileup.inframeIndels());
        assertEquals(2, (int) pileup.insertionCounts().get("TT"));
        assertEquals(18, (int) pileup.insertionCounts().get("CT"));
    }

    @Test
    public void testShortInsert() {
        final String line = "seq2\t156\tA\t11\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals("A", pileup.referenceBase());
        assertEquals(11, pileup.readCount());
        assertEquals(6, pileup.referenceCount());
        assertEquals(2, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(3, pileup.insertions());
        assertEquals(0, pileup.deletions());
        assertEquals(3, pileup.indels());
        assertEquals(3, (int) pileup.insertionCounts().get("AAG"));
    }

    @Test
    public void testLongInsert() {
        final String line = "seq2\t156\tA\t11\t.$......+2AG.+2AG.+12AGGGGGGGGGGGGG\t<975;:<<<<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(156, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(11, pileup.readCount());
        assertEquals(6, pileup.referenceCount());
        assertEquals(2, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(3, pileup.insertions());
        assertEquals(0, pileup.deletions());
        assertEquals(3, pileup.indels());
    }

    @Test
    public void testShortDelete() {
        final String line = "seq3\t200\tA\t20\t,,,,,..,.-4CACC.-4CACC....,.,,.^~.\t==<<<<<<<<<<<::<;2<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(200, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(20, pileup.readCount());
        assertEquals(18, pileup.referenceCount());
        assertEquals(0, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(0, pileup.cMismatchCount());
        assertEquals(0, pileup.insertions());
        assertEquals(2, pileup.deletions());
        assertEquals(2, pileup.indels());
        assertEquals(2, (int) pileup.deletionCounts().get("ACACC"));
    }

    @Test
    public void testLongDelete() {
        final String line = "seq3\t200\tA\t20\t,,,,,..,.-13CACCCCCCCCCCCC-4CACC....,.,,.^~.\t==<<<<<<<<<<<::<;2<<";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(200, pileup.position());
        assertEquals("A", pileup.referenceBase());
        assertEquals(20, pileup.readCount());
        assertEquals(18, pileup.referenceCount());
        assertEquals(0, pileup.gMismatchCount());
        assertEquals(0, pileup.aMismatchCount());
        assertEquals(0, pileup.tMismatchCount());
        assertEquals(1, pileup.cMismatchCount());
        assertEquals(0, pileup.insertions());
        assertEquals(2, pileup.deletions());
        assertEquals(2, pileup.indels());
    }
}
