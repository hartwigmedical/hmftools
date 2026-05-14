package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class LiftedMateInfoCacheTest
{
    @Test
    public void testGetPartnerReturnsOpposite()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();

        LiftedMateInfo r1 = LiftedMateInfo.mapped("1", 100, 199, "50M", false);
        LiftedMateInfo r2 = LiftedMateInfo.mapped("1", 400, 499, "50M", true);

        cache.recordPrimaryAlignment("read1", true, r1);
        cache.recordPrimaryAlignment("read1", false, r2);

        // looking up partner of R1 returns R2
        assertEquals(r2, cache.getPartnerMateInfo("read1", true));
        // looking up partner of R2 returns R1
        assertEquals(r1, cache.getPartnerMateInfo("read1", false));
    }

    @Test
    public void testGetPartnerReturnsNullWhenPartnerMissing()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        cache.recordPrimaryAlignment("read1", true, LiftedMateInfo.mapped("1", 100, 199, "50M", false));

        // R2 never recorded → partner-of-R1 is null
        assertNull(cache.getPartnerMateInfo("read1", true));
    }

    @Test
    public void testGetPartnerForUnknownReadIsNull()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        assertNull(cache.getPartnerMateInfo("missing", true));
    }

    @Test
    public void testUnmappedFactory()
    {
        LiftedMateInfo unmapped = LiftedMateInfo.unmapped();
        assertTrue(unmapped.Unmapped);
        assertNull(unmapped.Chromosome);
        assertEquals(0, unmapped.AlignmentStart);
        assertEquals(0, unmapped.AlignmentEnd);
    }
}
