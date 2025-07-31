package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class ProbeUtilsTest
{
    private static final int PROBE_START = 100;
    private static final int PROBE_CENTRE = 159;
    private static final int PROBE_END = 219;

    @Test
    public void testProbeRegionStartingAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionStartingAt(PROBE_START));
    }

    @Test
    public void testProbeRegionCenteredAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionCenteredAt(PROBE_CENTRE));
    }

    @Test
    public void testProbeRegionEndingAt()
    {
        assertEquals(new BaseRegion(PROBE_START, PROBE_END), probeRegionEndingAt(PROBE_END));
    }

    @Test
    public void testMinProbeStartContaining()
    {
        assertEquals(PROBE_START, ProbeUtils.minProbeStartContaining(PROBE_END));
    }

    @Test
    public void testMaxProbeEndContaining()
    {
        assertEquals(PROBE_END, ProbeUtils.maxProbeEndContaining(PROBE_START));
    }

    @Test
    public void testMinProbeStartOverlapping()
    {
        assertEquals(PROBE_START, minProbeStartOverlapping(new BaseRegion(PROBE_END, PROBE_END + 10)));
    }

    @Test
    public void testMaxProbeEndOverlapping()
    {
        assertEquals(PROBE_END, maxProbeEndOverlapping(new BaseRegion(PROBE_START - 10, PROBE_START)));
    }

    @Test
    public void testMinProbeStartWithoutGap()
    {
        assertEquals(PROBE_START, minProbeStartWithoutGap(new BaseRegion(PROBE_END + 1, PROBE_END + 10)));
    }

    @Test
    public void testMaxProbeEndWithoutGap()
    {
        assertEquals(PROBE_END, maxProbeEndWithoutGap(new BaseRegion(PROBE_START - 10, PROBE_START - 1)));
    }
}
