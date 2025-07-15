package com.hartwig.hmftools.geneutils.paneldesign;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ProbeUtilsTest
{
    private static final TargetRegion TARGET_REGION = new TargetRegion(
            new ChrBaseRegion("1", 100, 200),
            new TargetMetadata(TargetRegionType.CUSTOM, "extra"));
    private static final ProbeFactory PROBE_FACTORY = new ProbeFactory(TARGET_REGION);

    @Test
    public void testProbeCenteredAt()
    {
        CandidateProbe expected = new CandidateProbe(TARGET_REGION, new ChrBaseRegion("1", 90, 209));
        CandidateProbe actual = ProbeUtils.probeCenteredAt("1", 150, PROBE_FACTORY);
        assertEquals(expected, actual);
    }

    @Test
    public void testProbeStartingAt()
    {
        CandidateProbe expected = new CandidateProbe(TARGET_REGION, new ChrBaseRegion("1", 100, 219));
        CandidateProbe actual = ProbeUtils.probeStartingAt("1", 100, PROBE_FACTORY);
        assertEquals(expected, actual);
    }

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
