package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;

import org.junit.Test;

public class AlignmentScorerTest
{
    private static final String CHR = "chr1";

    // reference that returns all-'A' bases for any requested range.
    private static final RefSequenceSource ALL_A = (chromosome, posStart, posEnd) ->
    {
        byte[] bases = new byte[posEnd - posStart + 1];
        Arrays.fill(bases, (byte) 'A');
        return bases;
    };

    private static byte[] read(final String sequence)
    {
        return sequence.getBytes(StandardCharsets.US_ASCII);
    }

    private static int score(final String cigar, final byte[] readBases)
    {
        return AlignmentScorer.score(ALL_A, CHR, 1, cigar, readBases);
    }

    @Test
    public void testFullMatchScoresOnePerBase()
    {
        assertEquals(10, score("10M", read("AAAAAAAAAA")));
    }

    @Test
    public void testMismatchCostsFour()
    {
        // one C against an all-A reference: 9 matches (+9), 1 mismatch (-4) = 5.
        assertEquals(5, score("10M", read("AAAACAAAAA")));
    }

    @Test
    public void testDeletionIsAffineGap()
    {
        // 10 matched bases (+10) and a 2bp deletion (GAP_OPEN -6, GAP_EXTEND -1) = 3.
        assertEquals(3, score("5M2D5M", read("AAAAAAAAAA")));
    }

    @Test
    public void testIntronIsFree()
    {
        // the N contributes nothing; only the 10 matched bases score.
        assertEquals(10, score("5M100N5M", read("AAAAAAAAAA")));
    }

    @Test
    public void testSoftClipIsFree()
    {
        // 3 clipped bases contribute nothing; the 7 aligned bases score +7.
        assertEquals(7, score("3S7M", read("AAAAAAAAAA")));
    }

    @Test
    public void testNoRefSourceReturnsSentinel()
    {
        assertEquals(Integer.MIN_VALUE, AlignmentScorer.score(null, CHR, 1, "10M", read("AAAAAAAAAA")));
    }
}
