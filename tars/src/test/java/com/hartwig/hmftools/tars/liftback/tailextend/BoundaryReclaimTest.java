package com.hartwig.hmftools.tars.liftback.tailextend;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BoundaryReclaimTest
{
    @Test
    public void testMaxScoringPrefix()
    {
        // full match -> reclaim everything
        assertEquals(6, BoundaryReclaim.maxScoringPrefix(bases("ACGTAC"), bases("ACGTAC")));

        // leading mismatch (-4) with only 3 matches after (+3) never goes positive -> reclaim nothing
        assertEquals(0, BoundaryReclaim.maxScoringPrefix(bases("XCGT"), bases("ACGT")));

        // 5 matches then a trailing mismatch with no recovery -> cut at 5
        assertEquals(5, BoundaryReclaim.maxScoringPrefix(bases("ACGTAX"), bases("ACGTAC")));

        // one internal mismatch (-4) then enough matches to tie the peak -> reclaim through to the last match
        assertEquals(10, BoundaryReclaim.maxScoringPrefix(bases("AAAAAXAAAA"), bases("AAAAAAAAAA")));

        // mismatch after 5 matches, only 3 follow (+3 < +4 deficit) -> stays at the pre-mismatch peak
        assertEquals(5, BoundaryReclaim.maxScoringPrefix(bases("AAAAAXAAA"), bases("AAAAAAAAA")));

        // case-insensitive base comparison
        assertEquals(4, BoundaryReclaim.maxScoringPrefix(bases("acgt"), bases("ACGT")));
    }

    @Test
    public void testReversed()
    {
        assertArrayEquals(bases("TGCA"), BoundaryReclaim.reversed(bases("ACGT")));
    }
}
