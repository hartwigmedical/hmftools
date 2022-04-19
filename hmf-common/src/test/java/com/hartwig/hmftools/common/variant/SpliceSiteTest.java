package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.genome.region.Strand.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Strand.REVERSE;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SpliceSiteTest
{
    @Test
    public void testDonorSites()
    {
        assertEquals(-2, getDonorPosition(99, 100, FORWARD));
        assertEquals(-1, getDonorPosition(100, 100, FORWARD));
        assertEquals(1, getDonorPosition(101, 100, FORWARD));
        assertEquals(2, getDonorPosition(102, 100, FORWARD));
        assertEquals(5, getDonorPosition(105, 100, FORWARD));

        assertEquals(-2, getDonorPosition(101, 100, REVERSE));
        assertEquals(-1, getDonorPosition(100, 100, REVERSE));
        assertEquals(1, getDonorPosition(99, 100, REVERSE));
        assertEquals(2, getDonorPosition(98, 100, REVERSE));
        assertEquals(5, getDonorPosition(95, 100, REVERSE));
    }

    @Test
    public void testAcceptorSites()
    {
        assertEquals(3, getAcceptorPosition(97, 100, FORWARD));
        assertEquals(2, getAcceptorPosition(98, 100, FORWARD));
        assertEquals(1, getAcceptorPosition(99, 100, FORWARD));
        assertEquals(-1, getAcceptorPosition(100, 100, FORWARD));

        assertEquals(3, getAcceptorPosition(103, 100, REVERSE));
        assertEquals(2, getAcceptorPosition(102, 100, REVERSE));
        assertEquals(1, getAcceptorPosition(101, 100, REVERSE));
        assertEquals(-1, getAcceptorPosition(100, 100, REVERSE));
    }

}
