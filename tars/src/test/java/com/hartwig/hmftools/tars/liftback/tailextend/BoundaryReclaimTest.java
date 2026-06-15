package com.hartwig.hmftools.tars.liftback.tailextend;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.nio.charset.StandardCharsets;

import org.junit.Test;

public class BoundaryReclaimTest
{
    private static byte[] bases(final String text)
    {
        return text.getBytes(StandardCharsets.US_ASCII);
    }

    @Test
    public void testFullMatchReclaimsAll()
    {
        assertEquals(6, BoundaryReclaim.maxScoringPrefix(bases("ACGTAC"), bases("ACGTAC")));
    }

    @Test
    public void testFirstBaseMismatchNoRecoveryReclaimsNothing()
    {
        // leading mismatch (-4) with only 3 matches after (+3) never goes positive -> reclaim nothing
        assertEquals(0, BoundaryReclaim.maxScoringPrefix(bases("XCGT"), bases("ACGT")));
    }

    @Test
    public void testTrailingMismatchExcluded()
    {
        // 5 matches then a mismatch with no recovery -> cut at 5
        assertEquals(5, BoundaryReclaim.maxScoringPrefix(bases("ACGTAX"), bases("ACGTAC")));
    }

    @Test
    public void testInternalMismatchRecoveredIsReclaimed()
    {
        // one internal mismatch (-4) then enough matches to tie the peak -> reclaim through to the last match
        assertEquals(10, BoundaryReclaim.maxScoringPrefix(bases("AAAAAXAAAA"), bases("AAAAAAAAAA")));
    }

    @Test
    public void testInternalMismatchNotYetRecoveredCutsAtPeak()
    {
        // mismatch after 5 matches, only 3 matches follow (+3 < +4 deficit) -> stays at the pre-mismatch peak
        assertEquals(5, BoundaryReclaim.maxScoringPrefix(bases("AAAAAXAAA"), bases("AAAAAAAAA")));
    }

    @Test
    public void testCaseInsensitive()
    {
        assertEquals(4, BoundaryReclaim.maxScoringPrefix(bases("acgt"), bases("ACGT")));
    }

    @Test
    public void testReversed()
    {
        assertArrayEquals(bases("TGCA"), BoundaryReclaim.reversed(bases("ACGT")));
    }
}
