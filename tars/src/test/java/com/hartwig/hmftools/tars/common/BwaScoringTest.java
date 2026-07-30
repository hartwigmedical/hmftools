package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

// bwa-mem base scoring: MATCH +1, MISMATCH -4, case-insensitive.
public class BwaScoringTest
{
    @Test
    public void testScoreSumsMatchesAndMismatchesOverOverlap()
    {
        // 9 matches (+9) then 1 trailing mismatch (-4) over the 10-base overlap sum to +5.
        int result = BwaScoring.score(bases("A".repeat(10)), bases("A".repeat(9) + "C"));
        assertEquals(5, result);
    }

    @Test
    public void testScoreUsesOverlapLengthWhenInputsDiffer()
    {
        // min(10, 4) = 4-base overlap: the 4 leading matches score +4 and the 6 unpaired read bases are not scored.
        int result = BwaScoring.score(bases("A".repeat(10)), bases("AAAA"));
        assertEquals(4, result);
    }

    @Test
    public void testScoreTreatsBasesCaseInsensitively()
    {
        // lowercase read vs uppercase ref: basesEqual counts each as a match (+4).
        int result = BwaScoring.score(bases("aaaa"), bases("AAAA"));
        assertEquals(4, result);
    }

    @Test
    public void testMaxScoringPrefixReturnsFullLengthWhenAllMatch()
    {
        // every base matches, so the whole 10-base window is the highest-scoring prefix.
        int result = BwaScoring.maxScoringPrefix(bases("A".repeat(10)), bases("A".repeat(10)));
        assertEquals(10, result);
    }

    @Test
    public void testMaxScoringPrefixIsZeroWhenNoPrefixScoresPositive()
    {
        // every base mismatches, so the cumulative score never rises above 0 and nothing is extended.
        int result = BwaScoring.maxScoringPrefix(bases("C".repeat(5)), bases("A".repeat(5)));
        assertEquals(0, result);
    }

    @Test
    public void testMaxScoringPrefixStopsAtLongestPositivePrefix()
    {
        // 5 matches (peak +5) then 3 mismatches (dropping to +1, -3, -7): the prefix at the peak is the leading 5.
        int result = BwaScoring.maxScoringPrefix(bases("A".repeat(5) + "C".repeat(3)), bases("A".repeat(8)));
        assertEquals(5, result);
    }

    @Test
    public void testMaxScoringPrefixExtendsMatchesAfterInternalMismatch()
    {
        // 5 matches, 1 mismatch, 4 matches vs all-A: the >= tie climbs back to the +5 peak so the prefix extends to the full 10.
        int result = BwaScoring.maxScoringPrefix(bases("A".repeat(5) + "C" + "A".repeat(4)), bases("A".repeat(10)));
        assertEquals(10, result);
    }
}
