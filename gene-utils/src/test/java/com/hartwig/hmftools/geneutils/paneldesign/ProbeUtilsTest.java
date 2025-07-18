package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ProbeUtilsTest
{
    @Test
    public void testMinProbeStartContaining()
    {
        assertEquals(100, ProbeUtils.minProbeStartContaining(219));
    }

    @Test
    public void testMaxProbeEndContaining()
    {
        assertEquals(219, ProbeUtils.maxProbeEndContaining(100));
    }

    @Test
    public void testNextProbeStartPosition()
    {
        // Max overlap: 200, 199, 198, 197, 196, 195, 194, 193, 192, 191
        assertEquals(191, ProbeUtils.nextProbeStartPosition(200));
    }
}
