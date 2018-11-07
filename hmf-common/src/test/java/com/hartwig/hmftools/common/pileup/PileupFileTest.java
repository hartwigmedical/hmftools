package com.hartwig.hmftools.common.pileup;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import org.jetbrains.annotations.NotNull;
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
        assertEquals(25, pileup.insertCount());
        assertEquals(4, pileup.deleteCount());
        assertEquals(29, pileup.indelCount());
        assertEquals(1, pileup.inframeInserts().size());
        assertEquals(1, pileup.deleteCount("CTTT"));
        assertEquals(2, pileup.insertCount("TT"));
        assertEquals(18, pileup.insertCount("CT"));
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
        assertEquals(3, pileup.insertCount());
        assertEquals(0, pileup.deleteCount());
        assertEquals(3, pileup.indelCount());
        assertEquals(3, pileup.insertCount("AAG"));
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
        assertEquals(3, pileup.insertCount());
        assertEquals(0, pileup.deleteCount());
        assertEquals(3, pileup.indelCount());
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
        assertEquals(0, pileup.insertCount());
        assertEquals(2, pileup.deleteCount());
        assertEquals(2, pileup.indelCount());
        assertEquals(2, pileup.deleteCount("ACACC"));
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
        assertEquals(0, pileup.insertCount());
        assertEquals(2, pileup.deleteCount());
        assertEquals(2, pileup.indelCount());
    }

    @Test
    public void testQualScore() {
        assertEquals(0, PileupFile.qualityScore(0, "!"));
        assertEquals(2, PileupFile.qualityScore(0, "#"));
        assertEquals(3, PileupFile.qualityScore(0, "$"));

        final String line = "seq3\t200\tG\t12\t.,AaTtCc.+1Tt+1t.-1TT-1t^~\t!$#$##!!#$$#";
        final Pileup pileup = PileupFile.fromString(line);
        assertEquals(2, pileup.referenceCount());
        assertEquals(3, pileup.referenceScore());

        assertEquals(0, pileup.mismatchCount('G'));
        assertEquals(0, pileup.mismatchScore('G'));

        assertEquals(2, pileup.mismatchCount('A'));
        assertEquals(5, pileup.mismatchScore('A'));

        assertEquals(4, pileup.mismatchCount('T'));
        assertEquals(9, pileup.mismatchScore('T'));

        assertEquals(2, pileup.mismatchCount('C'));
        assertEquals(0, pileup.mismatchScore('C'));

        assertEquals(1, pileup.insertCount("GT"));
        assertEquals(2, pileup.insertScore("GT"));

        assertEquals(1, pileup.insertCount("TT"));
        assertEquals(3, pileup.insertScore("TT"));

        assertEquals(1, pileup.deleteCount("GT"));
        assertEquals(3, pileup.deleteScore("GT"));

        assertEquals(1, pileup.deleteCount("TT"));
        assertEquals(2, pileup.deleteScore("TT"));
    }

    @Test
    public void testQualScoreInterpreation() {
        final String qualities = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        for (int i = 0; i < qualities.length(); i++) {
            assertEquals(i, PileupFile.qualityScore(i, qualities));
        }
    }

    @NotNull
    public static Pileup create(@NotNull final String chromosome, final long position, int readCount, @NotNull final String ref,
            @NotNull final Map<String, Integer> counts, @NotNull final Map<String, Integer> score) {
        return ImmutablePileup.builder()
                .chromosome(chromosome)
                .position(position)
                .readCount(readCount)
                .referenceBase(ref)
                .countMap(counts)
                .scoreMap(score)
                .build();
    }

}
